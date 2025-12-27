#pragma once


class FramebufferRenderer {
public:

    explicit FramebufferRenderer(void* sdl_window, int max_frames_queue, bool vertical_flip = false);

    void push_frame(int* framebuffer, volatile bool* ready_flag) const;

    FramebufferRenderer(const FramebufferRenderer&) = delete;

    FramebufferRenderer operator=(const FramebufferRenderer&) = delete;

    ~FramebufferRenderer();

private:

    void* m_pimpl;

    friend FramebufferRenderer create_framebuffer_display(void* sdl_window);
};

