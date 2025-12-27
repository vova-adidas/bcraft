#pragma once

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>


#include <span>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <memory>
#include <glm/glm.hpp>
#include <bgl/bgl.hpp>

#include "../../concurrency/concurrent_queue.hpp"
#include "../../concurrency/parallel_for.hpp"

#include "../view.hpp"

#include "render_target.hpp"
#include "texture_rgba.hpp"
#include "warp.hpp"


inline bgl::ClipSpaceFrustum<bgl::ClipSpaceStyle::OpenGL> frustrum;

struct ClipVertex {
    float x, y, z, w, w_0, w_1;

    ClipVertex lerp(const ClipVertex& other, float t) const {
        return {
            x + (other.x - x) * t,
            y + (other.y - y) * t,
            z + (other.z - z) * t,
            w + (other.w - w) * t,
            w_0 + (other.w_0 - w_0) * t,
            w_1 + (other.w_1 - w_1) * t,
        };
    }
};

inline std::array<glm::vec2, 3> get_face_uv(const std::array<glm::vec3, 3>& v) {
    std::array<glm::vec2, 3> uv {};

    glm::vec3 e1{ v[1].x - v[0].x, v[1].y - v[0].y, v[1].z - v[0].z };
    glm::vec3 e2{ v[2].x - v[0].x, v[2].y - v[0].y, v[2].z - v[0].z };

    glm::vec3 n{
        e1.y * e2.z - e1.z * e2.y,
        e1.z * e2.x - e1.x * e2.z,
        e1.x * e2.y - e1.y * e2.x
    };

    float absX = std::fabs(n.x);
    float absY = std::fabs(n.y);
    float absZ = std::fabs(n.z);

    for (int i = 0;i < 3;i++) {
        if (absX >= absY && absX >= absZ)
            uv[i] = { v[i].z, v[i].y };
        else if (absY >= absX && absY >= absZ)
            uv[i] = { v[i].x, v[i].z };
        else
            uv[i] = { v[i].x, v[i].y };
    }

    return uv;
}

struct ChunkRenderUnit {
    std::vector<Vertex> vertices{};
    std::vector<int> indices{};
    glm::vec3 min{}, max{};
    int x{}, y{}, z{};
};

class SoftwareRenderer : public Renderer {
public:

    std::vector<std::unique_ptr<ChunkRenderUnit>> render_units {};
    std::vector<std::pair<float, ChunkRenderUnit*>> sorted_units_buffer {};

    ConcurrentQueue<std::move_only_function<void()>> run_render_thread {};

    // int* texture;
    // int texture_width{}, texture_height{};
    RenderTarget render_target = RenderTarget(0, 0);

    TextureRgba texture;
    int* tex_data;
    int tex_width, tex_height, channels;

    SoftwareRenderer() : texture(TextureRgba(nullptr, 0, 0)) {
        tex_data = reinterpret_cast<int*>(stbi_load(R"(texture.png)", &tex_width, &tex_height, &channels, 4));
        if (!tex_data)
            throw std::runtime_error("ТЕКСТУРА НЕ НАЙДЕНА");

        texture = TextureRgba(tex_data, tex_width, tex_height);
    }

    void* new_render_unit(
        std::shared_ptr<MeshData> data,
        int x,
        int y,
        int z,
        int local_min_x,
        int local_min_y,
        int local_min_z,
        int local_max_x,
        int local_max_y,
        int local_max_z) override {

        auto unit = std::make_unique<ChunkRenderUnit>();

        unit->x = x;
        unit->y = y;
        unit->z = z;
        unit->min = { local_min_x, local_min_y, local_min_z };
        unit->max = { local_max_x, local_max_y, local_max_z };
        unit->vertices = std::vector<Vertex>(data->vertices, data->vertices + data->num_vertices);
        unit->indices = std::vector<int>(data->indices, data->indices + data->num_indices);

        auto addr = unit.get();
        run_render_thread.produce([&, unit = std::move(unit)]() mutable {
            render_units.emplace_back(std::move(unit));
            sorted_units_buffer.emplace_back(0, nullptr);
        });

        return addr;

    }

    void delete_render_unit(void* unit_addr) override {
        run_render_thread.produce([&, unit_addr]() mutable {
            auto count = std::erase_if(render_units, [unit_addr](auto& it) {
                return it.get() == unit_addr;
            });
            
            for(std::size_t i = 0; i < count; ++i)
                sorted_units_buffer.pop_back();
        });
    }

    static std::pair<std::vector<glm::vec3>, std::vector<int>> aabb_mesh(const glm::vec3& min, const glm::vec3& max) {

        std::vector<glm::vec3> vertices{};
        std::vector<int> indices{};

        glm::vec3 v[8] = {
            {min.x, min.y, min.z},
            {max.x, min.y, min.z},
            {max.x, max.y, min.z},
            {min.x, max.y, min.z},
            {min.x, min.y, max.z},
            {max.x, min.y, max.z},
            {max.x, max.y, max.z},
            {min.x, max.y, max.z}
        };

        int faces[6][4] = {
            {0, 1, 2, 3},
            {5, 4, 7, 6},
            {4, 0, 3, 7},
            {1, 5, 6, 2},
            {4, 5, 1, 0},
            {3, 2, 6, 7}
        };

        for (auto & face : faces) {
            auto baseIndex = static_cast<int>(vertices.size());

            for (int i : face)
                vertices.push_back(v[i]);

            indices.insert(indices.end(), {
                baseIndex + 2, baseIndex + 1, baseIndex + 0,
                baseIndex + 3, baseIndex + 2, baseIndex + 0
            });
        }

        return { std::move(vertices), std::move(indices) };
    }

    bool depth_cull(ChunkRenderUnit* unit, const glm::mat4& pvm) {
        bgl::Viewport viewport(render_target.width(), render_target.height());

        auto min = unit->min;
        auto max = unit->max;

        auto [v_aabb, i_aabb] = aabb_mesh(min, max);
        for (std::size_t j = 0; j < std::size(i_aabb); j += 3) {
            auto aa = pvm * glm::vec4(v_aabb[i_aabb[j]], 1);
            auto ab = pvm * glm::vec4(v_aabb[i_aabb[j + 1]], 1);
            auto ac = pvm * glm::vec4(v_aabb[i_aabb[j + 2]], 1);

            auto edge = [](auto ax, auto ay, auto bx, auto by, auto x, auto y) {
                return (by - ay) * x + (ax - bx) * y + (bx * ay - ax * by);
                };

            ClipVertex vs[10];
            vs[0] = { aa.x, aa.y, aa.z, aa.w, 0.f, 0.f };
            vs[1] = { ab.x, ab.y, ab.z, ab.w, 0.f, 0.f };
            vs[2] = { ac.x, ac.y, ac.z, ac.w, 0.f, 0.f };

            int num = 3;
            frustrum.clip_polygon(vs, num);

            if (num < 3) continue;

            auto max_inv_w = std::numeric_limits<float>::lowest();
            for (int i = 0; i < num; ++i) {
                auto& v = vs[i];
                viewport.clip_space_to_screen_inv_w(v.x, v.y, v.w);
                if (v.w > max_inv_w)
                    max_inv_w = v.w;
            }

            for (int i = 1; i < num - 1; ++i) {
                auto a = vs[0];
                auto b = vs[i];
                auto c = vs[i + 1];

                if (edge(a.x, a.y, b.x, b.y, c.x, c.y) <= 0)
                    continue;

                bgl::TriangleScanlineIterator scanline(1.0f * a.x / TILE_SIZE, 1.0f * a.y / TILE_SIZE, 1.0f * b.x / TILE_SIZE, 1.0f * b.y / TILE_SIZE, 1.0f * c.x / TILE_SIZE, 1.0f * c.y / TILE_SIZE);

                while (scanline.next()) {

                    int y = scanline.line();
                    for (int x = scanline.begin(); x < scanline.end(); ++x)
                        if (max_inv_w > render_target.read_depth_mipmap(x, y))
                            return false;
                }
            }
        }
        return true;
    }


    static float nearest_distance(const glm::vec3& min_view, const glm::vec3& max_view) {
    
        //auto center = (min_view + max_view) / 2.0f;

        //return glm::dot(center, center);

        glm::vec3 nearest;
        nearest.x = 0 < min_view.x ? min_view.x : (0 > max_view.x ? max_view.x : 0);
        nearest.y = 0 < min_view.y ? min_view.y : (0 > max_view.y ? max_view.y : 0);
        nearest.z = 0 < min_view.z ? min_view.z : (0 > max_view.z ? max_view.z : 0);
        return -nearest.z; 
    }

    int frustum_cull_and_sort_units(const glm::mat4& projection, const glm::mat4& view) {
        auto projection_view = projection * view;

        int active_units = 0;

        for (auto& unit : render_units) {

            auto min = unit->min;
            auto max = unit->max;

            if (min == max)
                continue;

            auto model = glm::translate(glm::identity<glm::mat4>(), glm::vec3(unit->x, unit->y, unit->z));

            glm::vec3 min_view = view * model * glm::vec4(min, 1);
            glm::vec3 max_view = view * model * glm::vec4(max, 1);

            float distance_to_camera = nearest_distance(min_view, max_view);

            auto pvm = projection_view * model;

            bgl::ClipSpaceFrustum<bgl::ClipSpaceStyle::OpenGL> frustrum;

            glm::vec4 corners[8] = {
                pvm * glm::vec4(min.x, min.y, min.z, 1.0f),
                pvm * glm::vec4(max.x, min.y, min.z, 1.0f),
                pvm * glm::vec4(max.x, max.y, min.z, 1.0f),
                pvm * glm::vec4(min.x, max.y, min.z, 1.0f),
                pvm * glm::vec4(min.x, min.y, max.z, 1.0f),
                pvm * glm::vec4(max.x, min.y, max.z, 1.0f),
                pvm * glm::vec4(max.x, max.y, max.z, 1.0f),
                pvm * glm::vec4(min.x, max.y, max.z, 1.0f)
            };

            if (frustrum.outside(corners, 8))
                continue;

            sorted_units_buffer[active_units++] = { distance_to_camera , unit.get() };
        }

        std::sort(std::begin(sorted_units_buffer), std::begin(sorted_units_buffer) + active_units, [](auto& lhs, auto& rhs) { return lhs.first < rhs.first; });

        return active_units;
    }

    struct PrimitiveBuffer {
        constexpr static int Size = 2'000'000;

        alignas(Warp::ALIGNMENT) float light_0[Size];
        alignas(Warp::ALIGNMENT) float light_1[Size];
        alignas(Warp::ALIGNMENT) float light_2[Size];
        alignas(Warp::ALIGNMENT) float u_0[Size];
        alignas(Warp::ALIGNMENT) float v_0[Size];
        alignas(Warp::ALIGNMENT) float u_1[Size];
        alignas(Warp::ALIGNMENT) float v_1[Size];
        alignas(Warp::ALIGNMENT) float u_2[Size];
        alignas(Warp::ALIGNMENT) float v_2[Size];
        alignas(Warp::ALIGNMENT) int tex[Size];
        alignas(Warp::ALIGNMENT) float fog_0[Size];
        alignas(Warp::ALIGNMENT) float fog_1[Size];
        alignas(Warp::ALIGNMENT) float fog_2[Size];
    };

    PrimitiveBuffer& raster_to_vis_buffer(auto units_begin, auto units_end, RenderTarget& render_target, float render_distance, const glm::mat4& projection, const glm::mat4& view) {

        static PrimitiveBuffer primitive_buffer;

        std::atomic<int> prim_index;
        std::atomic<int> total_triangles;
        std::atomic<int> culled;
        std::atomic<int> current = -1;
        parallel_foreach(units_begin, units_end, [&](auto pair) {
            auto th = current.fetch_add(1);
            ChunkRenderUnit* unit = pair.second;

            auto model = glm::translate(glm::identity<glm::mat4>(), glm::vec3(unit->x, unit->y, unit->z)) *
                glm::scale(glm::identity<glm::mat4>(), glm::vec3(1.f, 1.f, 1.f));

            bgl::Viewport viewport(render_target.width(), render_target.height());

            auto pvm = projection * view *
                glm::translate(glm::identity<glm::mat4>(), glm::vec3(unit->x, unit->y, unit->z)) *
                glm::scale(glm::identity<glm::mat4>(), glm::vec3(1.f, 1.f, 1.f));

            if (depth_cull(unit, pvm)) {
                culled++;
                return;
            }

            int tris_count = std::size(unit->indices) / 3;
            total_triangles += tris_count;
            for (std::size_t i = 0; i < std::size(unit->indices); i += 3) {
                auto a = unit->vertices[unit->indices[i]];
                auto b = unit->vertices[unit->indices[i + 1]];
                auto c = unit->vertices[unit->indices[i + 2]];

                auto a_view = view * model * glm::vec4(a.pos, 1);
                auto b_view = view * model * glm::vec4(b.pos, 1);
                auto c_view = view * model * glm::vec4(c.pos, 1);

                auto apos = pvm * glm::vec4(a.pos, 1);
                auto bpos = pvm * glm::vec4(b.pos, 1);
                auto cpos = pvm * glm::vec4(c.pos, 1);

                float al = std::min(1.0f, (1.0f - a.ao / 4.0f) * pow(a.lighting / 15.0f, 1.4f));
                float bl = std::min(1.0f, (1.0f - b.ao / 4.0f) * pow(b.lighting / 15.0f, 1.4f));
                float cl = std::min(1.0f, (1.0f - c.ao / 4.0f) * pow(c.lighting / 15.0f, 1.4f));

                ClipVertex vertices[10];
                vertices[0] = { apos.x, apos.y, apos.z, apos.w, 1.f, 0.f };
                vertices[1] = { bpos.x, bpos.y, bpos.z, bpos.w, 0.f, 1.f };
                vertices[2] = { cpos.x, cpos.y, cpos.z, cpos.w, 0.f, 0.f };

                int num_vertices = 3;
                frustrum.clip_polygon(vertices, num_vertices);

                if (num_vertices < 3)
                    continue;

                for (int i = 0; i < num_vertices; ++i) {
                    auto& v = vertices[i];
                    viewport.clip_space_to_screen_inv_w(v.x, v.y, v.w);
                }

                int id = prim_index.fetch_add(1, std::memory_order_relaxed);
                if (id >= primitive_buffer.Size)
                    return;

                primitive_buffer.light_0[id] = al;
                primitive_buffer.light_1[id] = bl;
                primitive_buffer.light_2[id] = cl;

                auto uv = get_face_uv({ a.pos, b.pos, c.pos });

                primitive_buffer.u_0[id] = uv[0].x;
                primitive_buffer.v_0[id] = uv[0].y;
                primitive_buffer.u_1[id] = uv[1].x;
                primitive_buffer.v_1[id] = uv[1].y;
                primitive_buffer.u_2[id] = uv[2].x;
                primitive_buffer.v_2[id] = uv[2].y;

                primitive_buffer.tex[id] = a.texture;

                primitive_buffer.fog_0[id] = 1.0f - std::clamp(-a_view.z / ((render_distance - 2) * Chunk_width - 8), 0.0f, 1.0f);
                primitive_buffer.fog_1[id] = 1.0f - std::clamp(-b_view.z / ((render_distance - 2) * Chunk_width - 8), 0.0f, 1.0f);
                primitive_buffer.fog_2[id] = 1.0f - std::clamp(-c_view.z / ((render_distance - 2) * Chunk_width - 8), 0.0f, 1.0f);

                for (int i = 1; i < num_vertices - 1; ++i) {

                    auto v0 = vertices[0];
                    auto v1 = vertices[i];
                    auto v2 = vertices[i + 1];

                    constexpr int Num_attrs = 4;
                    constexpr int Weight_0 = 0;
                    constexpr int Weight_1 = 1;
                    constexpr int Weight_2 = 2;
                    constexpr int Inv_w = 3;

                    float a_attrs[Num_attrs]{
                        v0.w_0 * v0.w,
                        v0.w_1 * v0.w,
                        (1.0f - v0.w_1 - v0.w_0) * v0.w,
                        v0.w
                    };

                    float b_attrs[Num_attrs]{
                        v1.w_0 * v1.w,
                        v1.w_1 * v1.w,
                        (1.0f - v1.w_1 - v1.w_0) * v1.w,
                        v1.w
                    };

                    float c_attrs[Num_attrs]{
                        v2.w_0 * v2.w,
                        v2.w_1 * v2.w,
                        (1.0f - v2.w_1 - v2.w_0) * v2.w,
                        v2.w
                    };

                    float dmin = std::numeric_limits<float>::max();

                    bgl::triangle<Warp, TILE_SIZE>(
                        v0.x, v0.y,
                        v1.x, v1.y,
                        v2.x, v2.y,
                        a_attrs,
                        b_attrs,
                        c_attrs,
                        [&](int tx, int ty, int column, int row, Warp::IntMask mask, Warp::Float* attrs) {
                            int x = column + tx * TILE_SIZE;
                            int y = row + ty * TILE_SIZE;

                            auto prev_inv_w = render_target.read_depth(x, y);

                            auto depth_mask = attrs[Inv_w] > prev_inv_w;
                            mask &= xsimd::batch_bool_cast<int, float>(depth_mask);

                            auto fmask = xsimd::batch_bool_cast<float, int>(mask);

                            auto new_inv_w = xsimd::select(fmask, attrs[Inv_w], prev_inv_w);
                            render_target.write_depth(x, y, new_inv_w);

                            dmin = std::min(dmin, xsimd::reduce_min(new_inv_w));

                            render_target.write_weight_0(x, y, xsimd::select(fmask, attrs[Weight_0], render_target.read_weight_0(x, y)));
                            render_target.write_weight_1(x, y, xsimd::select(fmask, attrs[Weight_1], render_target.read_weight_1(x, y)));

                            render_target.write_primitive(x, y, xsimd::select(mask, Warp::Int(id), render_target.read_primitive(x, y)));

                        },
                        [&](int tx, int ty) {
                            dmin = std::numeric_limits<float>::max();
                            render_target.lock_tile(tx, ty);
                        },
                        [&](int tx, int ty) {
                            render_target.write_depth_mipmap(tx, ty, dmin);
                            render_target.unlock_tile(tx, ty);
                        }
                    );
                }
            }
        });

        return primitive_buffer;
    }

    void shade(int* framebuffer, const PrimitiveBuffer& primitive_buffer, const RenderTarget& render_target) {
         static std::vector<int> indices;
        indices.clear();
        indices.reserve(render_target.height());
        for (int i = 0; i < render_target.height(); ++i)
            indices.emplace_back(i);

        parallel_foreach(std::begin(indices), std::end(indices), [&](auto y) {
            for (int x = 0; x < render_target.width(); x += Warp::SIZE) {

                auto pid = render_target.read_primitive(x, y);

                auto mask = pid != -1;

                auto idx = xsimd::select(mask, pid, Warp::scalar(0));

                auto l0 = Warp::Float::gather(primitive_buffer.light_0, idx);
                auto l1 = Warp::Float::gather(primitive_buffer.light_1, idx);
                auto l2 = Warp::Float::gather(primitive_buffer.light_2, idx);
                auto u0 = Warp::Float::gather(primitive_buffer.u_0, idx);
                auto u1 = Warp::Float::gather(primitive_buffer.u_1, idx);
                auto u2 = Warp::Float::gather(primitive_buffer.u_2, idx);
                auto v0 = Warp::Float::gather(primitive_buffer.v_0, idx);
                auto v1 = Warp::Float::gather(primitive_buffer.v_1, idx);
                auto v2 = Warp::Float::gather(primitive_buffer.v_2, idx);
                auto f0 = Warp::Float::gather(primitive_buffer.fog_0, idx);
                auto f1 = Warp::Float::gather(primitive_buffer.fog_1, idx);
                auto f2 = Warp::Float::gather(primitive_buffer.fog_2, idx);
                auto tex_id = Warp::Int::gather(primitive_buffer.tex, idx);

                auto inv_w = render_target.read_depth(x, y);
                auto w0 = render_target.read_weight_0(x, y) / inv_w;
                auto w1 = render_target.read_weight_1(x, y) / inv_w;
                auto w2 = 1.0f - w0 - w1;
                auto l = l0 * w0 + l1 * w1 + l2 * w2;
                auto f = f0 * w0 + f1 * w1 + f2 * w2;

                auto u = u0 * w0 + u1 * w1 + u2 * w2;
                auto v = v0 * w0 + v1 * w1 + v2 * w2;

                float step = 16.0f / 255.f;

                u = u - xsimd::floor(u);
                v = Warp::scalar(1.0f) - (v - xsimd::floor(v));
                u *= step;
                v *= step;
                u += Warp::scalar(step) * xsimd::to_float(tex_id);

                Warp::Float r, g, b, a;
                texture.sample<Warp>(u, v, r, g, b, a);

                auto prev_color_r = Warp::scalar(1.0f);
                auto prev_color_g = Warp::scalar(1.0f);
                auto prev_color_b = Warp::scalar(1.0f);

                r *= l;
                g *= l;
                b *= l;


                auto fmask = xsimd::batch_bool_cast<float, int>(mask);
                auto final_r = xsimd::select(fmask, xsimd::clip(prev_color_g + (r - prev_color_g) * f, Warp::scalar(0.0f), Warp::scalar(1.0f)), prev_color_r);
                auto final_g = xsimd::select(fmask, xsimd::clip(prev_color_g + (g - prev_color_g) * f, Warp::scalar(0.0f), Warp::scalar(1.0f)), prev_color_g);
                auto final_b = xsimd::select(fmask, xsimd::clip(prev_color_b + (b - prev_color_b) * f, Warp::scalar(0.0f), Warp::scalar(1.0f)), prev_color_b);

                auto red = xsimd::to_int(final_r * Warp::scalar(255.f) + Warp::scalar(.5f));
                auto green = xsimd::to_int(final_g * Warp::scalar(255.f) + Warp::scalar(.5f));
                auto blue = xsimd::to_int(final_b * Warp::scalar(255.f) + Warp::scalar(.5f));

                auto red_shifted = red;
                auto green_shifted = green << 8;
                auto blue_shifted = blue << 16;
                auto alpha = Warp::Int(0xFF000000);

                auto rgba = red_shifted | green_shifted | blue_shifted | alpha;

                rgba.store_aligned(framebuffer + x + y * render_target.width());
            }
        });
    }

    void render(int* framebuffer, int fb_width, int fb_height, float render_distance, const glm::mat4& projection, const glm::mat4& view) {
        auto frame_begin = std::chrono::steady_clock::now();

        std::move_only_function<void()> task;
        while (run_render_thread.consume(task))
            task();

        render_target.resize(fb_width, fb_height);
        render_target.clear();

        int active_units = frustum_cull_and_sort_units(projection, view);

        auto begin = std::chrono::steady_clock::now();

        auto& primitive_buffer = raster_to_vis_buffer(std::begin(sorted_units_buffer), std::begin(sorted_units_buffer) + active_units, render_target, render_distance, projection, view);

        auto end = std::chrono::steady_clock::now();

        auto raster_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9;

        begin = std::chrono::steady_clock::now();

        shade(framebuffer, primitive_buffer, render_target);

        end = std::chrono::steady_clock::now();

        auto shade_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count()  / 1e9;

        auto frame_end = std::chrono::steady_clock::now();

        auto frame_time = std::chrono::duration_cast<std::chrono::nanoseconds>(frame_end - frame_begin).count()  / 1e9;

        std::cout << std::fixed << std::setprecision(5)
                << "[Raster: "  << std::setw(8) << std::left << raster_time << " sec.], "
                << "[Shading: " << std::setw(8) << std::left << shade_time  << " sec.], "
                << "[Fps: "     << std::setw(8) << std::left << (1.0 / frame_time) << "]"
                << std::endl;

    }
};