#pragma once

#include "warp.hpp"

class TextureRgba {
public:
    TextureRgba(int* data, int width, int height) : m_data(data), m_width(width), m_height(height) {

    }

    template<typename TWarp>
    void sample(
        TWarp::Float u, TWarp::Float v,
        TWarp::Float& out_r,
        TWarp::Float& out_g,
        TWarp::Float& out_b,
        TWarp::Float& out_a) {
        u = xsimd::clip(u, Warp::scalar(0.0f), Warp::scalar(1.0f));
        v = xsimd::clip(v, Warp::scalar(0.0f), Warp::scalar(1.0f));

        auto texture_x = xsimd::to_int(u * Warp::scalar((m_width - 1.0f)));
        auto texture_y = xsimd::to_int(v * Warp::scalar((m_height - 1.0f)));
        auto texture_idx = texture_x + texture_y * Warp::scalar(m_width);

        auto color = Warp::Int::gather(m_data, texture_idx);
        out_r = xsimd::to_float(color & 0x000000FF) / Warp::scalar(255.0f);
        out_g = xsimd::to_float(color >> 8 & 0x000000FF) / Warp::scalar(255.0f);
        out_b = xsimd::to_float(color >> 16 & 0x000000FF) / Warp::scalar(255.0f);
        out_a = xsimd::to_float(color >> 24 & 0x000000FF) / Warp::scalar(255.0f);
    }

    int width() const { return m_width; }

    int height() const { return m_height; }
private:
    int* m_data;
    int m_width, m_height;
};