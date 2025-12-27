#include "FramebufferRenderer.hpp"
#include <SDL.h>
#include <glad/glad.h>
#include <atomic>
#include <mutex>
#include <queue>
#include <thread>
#include <sstream>
#include <stdexcept>

struct FrameHandle {
    int* buffer;
    volatile bool* render_flag;
};

class FrameQueue {
public:
    explicit FrameQueue(std::size_t max_size) : max_size_(max_size) {}

    // Производитель
    void produce(FrameHandle frame) {
        std::unique_lock lock(mutex_);
        // Ждем, если очередь полна
        cv_not_full_.wait(lock, [this] { return queue_.size() < max_size_ || stop_; });
        if (stop_) return;

        queue_.push(frame);
        cv_not_empty_.notify_one();
    }

    // Потребитель
    bool consume(FrameHandle& frame) {
        std::unique_lock lock(mutex_);
        cv_not_empty_.wait(lock, [this]{ return !queue_.empty() || stop_; });
        if (queue_.empty())
            return false; // очередь пуста и остановлена

        frame = queue_.front();
        queue_.pop();
        cv_not_full_.notify_one();
        return true;
    }

    std::size_t size() {
        std::lock_guard lock(mutex_);
        return queue_.size();
    }

    void stop() {
        {
            std::lock_guard lock(mutex_);
            stop_ = true;
        }
        cv_not_empty_.notify_all();
        cv_not_full_.notify_all();
    }

private:
    std::queue<FrameHandle> queue_;
    std::mutex mutex_;
    std::condition_variable cv_not_empty_;
    std::condition_variable cv_not_full_;
    const std::size_t max_size_;
    bool stop_ = false;
};


class FramebufferDisplayImpl {
public:
    explicit FramebufferDisplayImpl(SDL_Window* window, int frames_queue_size, bool vertical_flip) :  m_queue(frames_queue_size), m_window(window), m_vertical_flip(vertical_flip) {

        auto gl_context = SDL_GL_CreateContext(window);
        m_gl_context = gl_context;
        if (!gl_context)
            throw std::runtime_error("Пошел нахуй + glContext error");
        SDL_GL_MakeCurrent(m_window, NULL);

        SDL_GetWindowSize(window, &m_width, &m_height);

        m_thread = std::thread([&]() {
            render();
        });

        std::unique_lock lock(m_mutex);
        m_init_wait.wait(lock, [this] { return m_initialized.load(); });
        if (m_failed)
            throw std::runtime_error(std::string(m_error_message));
    }

    ~FramebufferDisplayImpl() {
        m_queue.stop();
        m_thread.join();
        SDL_GL_DeleteContext(m_gl_context);
    }

    void push_frame(int* framebuffer, volatile bool* rendered) {
        m_queue.produce({framebuffer, rendered});
    }

    std::size_t num_frames() {
        return m_queue.size();
    }
private:

    FrameQueue m_queue;

    SDL_Window* m_window;
    SDL_GLContext m_gl_context {};
    std::thread m_thread;
    GLuint m_pbos[2] {};
    GLuint m_texture {};

    std::condition_variable m_init_wait {};
    std::mutex m_mutex {};

    int m_index {}, m_width {}, m_height {};

    bool m_vertical_flip;
    std::atomic<bool> m_failed {};
    std::atomic<bool> m_initialized = false;
    char m_error_message[100] {};

    bool init() {
        SDL_GL_MakeCurrent(m_window, m_gl_context);
        std::stringstream serr;

        auto set_error = [&] {
            auto string = serr.str();
            int i = 0;
            for (i = 0; i < string.length(); ++i) {
                if (i >= 99)
                    break;
                m_error_message[i] = string[i];
            }
            m_error_message[i] = '\0';
        };

        if (!gladLoadGLLoader((GLADloadproc)SDL_GL_GetProcAddress)) {
            serr << "Failed to initialize GLAD";
            set_error();
            return false;
        }

        SDL_GL_SetSwapInterval(1);

        glGenTextures(1, &m_texture);
        glBindTexture(GL_TEXTURE_2D, m_texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, m_width, m_height, 0,
           GL_RGBA, GL_UNSIGNED_BYTE, nullptr);

        glGenBuffers(2, m_pbos);
        for (int i = 0; i < 2; i++) {
            glBindBuffer(GL_PIXEL_UNPACK_BUFFER, m_pbos[i]);
            glBufferData(GL_PIXEL_UNPACK_BUFFER, m_width * m_height * 4, nullptr, GL_STREAM_DRAW);
        }
        glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);

        m_index = 0;
        return true;
    }

    void deinit() {
        glDeleteBuffers(2, m_pbos);
        glDeleteTextures(1, &m_texture);
    }

    void render() {
        m_failed = !init();

        {
            std::lock_guard lock(m_mutex);
            m_initialized = true;
        }
        m_init_wait.notify_one();

        if (m_failed)
            return;

        const auto drawQuad = [](bool flip) {
            if (flip) {
                glTexCoord2f(0,0); glVertex2f(-1,-1);
                glTexCoord2f(1,0); glVertex2f(1,-1);
                glTexCoord2f(1,1); glVertex2f(1,1);
                glTexCoord2f(0,1); glVertex2f(-1,1);
            }
            else {
                glTexCoord2f(0,1); glVertex2f(-1,-1);
                glTexCoord2f(1,1); glVertex2f(1,-1);
                glTexCoord2f(1,0); glVertex2f(1,1);
                glTexCoord2f(0,0); glVertex2f(-1,1);
            }
        };

        FrameHandle framebuffer{};
        if (m_queue.consume(framebuffer)) {
            SDL_GL_MakeCurrent(m_window, m_gl_context);
            // Жёсткая загрузка первого кадра
            glBindBuffer(GL_PIXEL_UNPACK_BUFFER, m_pbos[m_index]);
            glBufferData(GL_PIXEL_UNPACK_BUFFER, m_width * m_height * 4, nullptr, GL_STREAM_DRAW);
            void* ptr = glMapBuffer(GL_PIXEL_UNPACK_BUFFER, GL_WRITE_ONLY);
            if (ptr) {
                memcpy(ptr, framebuffer.buffer, m_width * m_height * 4);
                glUnmapBuffer(GL_PIXEL_UNPACK_BUFFER);
            }

            glBindTexture(GL_TEXTURE_2D, m_texture);
            glBindBuffer(GL_PIXEL_UNPACK_BUFFER, m_pbos[m_index]);
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, m_width, m_height,
                            GL_RGBA, GL_UNSIGNED_BYTE, nullptr);

            // Отрисовка
            glClear(GL_COLOR_BUFFER_BIT);
            glEnable(GL_TEXTURE_2D);
            glBegin(GL_QUADS);

            drawQuad(m_vertical_flip);

            glEnd();
            glDisable(GL_TEXTURE_2D);

            SDL_GL_SwapWindow(m_window);

            // Переходим к следующему индексу для двойной буферизации
            m_index = (m_index + 1) % 2;
            *framebuffer.render_flag = true;
        }

        while (m_queue.consume(framebuffer)) {
            SDL_GL_MakeCurrent(m_window, m_gl_context);
            int nextIndex = (m_index + 1) % 2;

            // Записываем данные в следующий PBO
            glBindBuffer(GL_PIXEL_UNPACK_BUFFER, m_pbos[nextIndex]);
            glBufferData(GL_PIXEL_UNPACK_BUFFER, m_width * m_height * 4, nullptr, GL_STREAM_DRAW);
            void* ptr = glMapBuffer(GL_PIXEL_UNPACK_BUFFER, GL_WRITE_ONLY);
            if (ptr) {
                memcpy(ptr, framebuffer.buffer, m_width * m_height * 4);
                glUnmapBuffer(GL_PIXEL_UNPACK_BUFFER);
            }

            // Используем предыдущий буфер для TexSubImage2D
            glBindTexture(GL_TEXTURE_2D, m_texture);
            glBindBuffer(GL_PIXEL_UNPACK_BUFFER, m_pbos[m_index]);
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, m_width, m_height,
                            GL_RGBA, GL_UNSIGNED_BYTE, nullptr);

            // Отрисовка
            glClear(GL_COLOR_BUFFER_BIT);
            glEnable(GL_TEXTURE_2D);
            glBegin(GL_QUADS);
            drawQuad(m_vertical_flip);
            glEnd();
            glDisable(GL_TEXTURE_2D);

            SDL_GL_SwapWindow(m_window);

            m_index = nextIndex;

            *framebuffer.render_flag = true;
        }

        deinit();
    }

};

FramebufferRenderer::FramebufferRenderer(void* sdl_window, int max_frames_queue, bool vertical_flip) : m_pimpl(new FramebufferDisplayImpl(static_cast<SDL_Window*>(sdl_window), max_frames_queue, vertical_flip)){

}


FramebufferRenderer::~FramebufferRenderer() {
    delete static_cast<FramebufferDisplayImpl*>(m_pimpl);
}

void FramebufferRenderer::push_frame(int *framebuffer, volatile bool* ready_flag) const {
    static_cast<FramebufferDisplayImpl*>(m_pimpl)->push_frame(framebuffer, ready_flag);
}
