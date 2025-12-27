#include <bgl/bgl.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <execution>
#include <SDL.h>
#include "core/world.hpp"
#include "display/FramebufferRenderer.hpp"
#include "rendering/software/warp.hpp"
#include "rendering/view.hpp"
#include "rendering/software/software_renderer.hpp"

constexpr int WINDOW_WIDTH = 1920;
constexpr int WINDOW_HEIGHT = 1080;
constexpr int RENDER_DISTANCE = 15;

class Camera {
public:
   glm::vec3 position{ 0.f, 0.f, 0.f };
   glm::vec3 target{ 0.f, 0.f, -1.f };
   glm::vec3 up{ 0.f, 1.f, 0.f };

   Camera(float fov_y_rad, float aspect, float near, float far) : m_projection(glm::perspective(fov_y_rad, aspect, near, far)) {}

   glm::mat4x4 projection_view() const { return m_projection * view(); }

   glm::mat4x4 view() const { return glm::lookAt(position, target, up); }

   const glm::mat4x4& projection() const { return m_projection; }

private:

   glm::mat4x4 m_projection;
};

class Client : public Actor {
public:
    FramebufferRenderer* display;
    std::shared_ptr<SoftwareRenderer> renderer;
    View view;
    float render_distance;

    Client(FramebufferRenderer* display, float render_distance) : renderer(std::make_shared<SoftwareRenderer>()), display(display), view(renderer), render_distance(render_distance) {
    }

    void receive_chunk(std::shared_ptr<Chunk> chunk) override {
        view.push_chunk(chunk);

    }

    void update() {
        view.update();
    }

    int ring {};

    void render(const glm::mat4& projection, const glm::mat4& m_view) const {

        alignas(Warp::ALIGNMENT) static int framebuffers[4][WINDOW_WIDTH * WINDOW_HEIGHT];
        static volatile bool available[4] { true, true, true, true };

        int index = -1;
        while (index == -1)
            for (int i = 0; i < 4; ++i)
                if (available[i]) {
                    index = i;
                    break;
                }

        auto framebuffer = framebuffers[index];

        renderer->render(framebuffer, WINDOW_WIDTH, WINDOW_HEIGHT, render_distance, projection, m_view);
        available[index] = false;
        display->push_frame(framebuffer, available + index);

    }
};


int main(int, char**) {
	using namespace std;

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "SDL Init error: " << SDL_GetError() << std::endl;
        return -1;
    }

    SDL_Window* window = SDL_CreateWindow("SDL + Double PBO",
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        WINDOW_WIDTH,
        WINDOW_HEIGHT,
        SDL_WINDOW_OPENGL);

    if (!window) {
        std::cerr << "SDL CreateWindow error: " << SDL_GetError() << std::endl;
        return -1;
    }

    constexpr int frameQueueSize = 4;
    bool flip = true;
    FramebufferRenderer renderer(window, frameQueueSize, flip);

    Camera cam(glm::radians(75.f), 1.0f * WINDOW_WIDTH / WINDOW_HEIGHT, 0.01f, 1000.0f);
    cam.position = { 0, 0, 0 };
    auto camForward = glm::normalize(cam.target - cam.position);
    glm::vec3 camUp = { 0,1,0 };

    float yaw = 0.0f;
    float pitch = 0.0f;
    float moveSpeed = 1.0f;
    float mouseSensitivity = 0.002f;

    World world(RENDER_DISTANCE);

    Client client(&renderer, world.load_distance());
    world.add_actor(&client);

    SDL_Event event;
    bool running = true;

    while (running) {

        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) running = false;
            if (event.type == SDL_MOUSEMOTION && (SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(SDL_BUTTON_LEFT))) {
                yaw += event.motion.xrel * mouseSensitivity;
                pitch -= event.motion.yrel * mouseSensitivity;

                if (pitch > M_PI / 2 - 0.01f) pitch = M_PI / 2 - 0.01f;
                if (pitch < -M_PI / 2 + 0.01f) pitch = -M_PI / 2 + 0.01f;

                camForward.x = glm::cos(pitch) * glm::sin(yaw);
                camForward.y = glm::sin(pitch);
                camForward.z = -glm::cos(pitch) * glm::cos(yaw);
                camForward = glm::normalize(camForward);

                cam.target = cam.position + camForward;
            }

            if (event.type == SDL_MOUSEWHEEL) {
                cam.position += camForward * (float)event.wheel.y * moveSpeed;
                cam.target = cam.position + camForward;

                
            }
        }

        client.pos.x = cam.position.x;
        client.pos.y = cam.position.y;
        client.pos.z = cam.position.z;

        world.tick();
        client.update();

        client.render(cam.projection(), cam.view());
    

    }


    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
