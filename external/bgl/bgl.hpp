#pragma once

#include "clipping.hpp"
#include "raster.hpp"
#include "viewport.hpp"

namespace bgl {


    template<typename TWarp>
    struct Config {
        using Warp = TWarp;

        ClipSpaceStyle clip_space_style;
        int tile_size;
        bool perspective_correction;
        bool use_reciprocal_for_perspective_correction;

        inline constexpr Config() {};

        inline constexpr Config(int tile_size, ClipSpaceStyle clip_space_style, bool perspective_correction, bool use_reciprocal_for_perspective_correction)
            : clip_space_style(clip_space_style), tile_size(tile_size), perspective_correction(perspective_correction), use_reciprocal_for_perspective_correction(use_reciprocal_for_perspective_correction) {}

        inline constexpr Config(int tile_size, ClipSpaceStyle clip_space_style, bool perspective_correction)
            : clip_space_style(clip_space_style), tile_size(tile_size), perspective_correction(perspective_correction), use_reciprocal_for_perspective_correction(false) {}

        inline constexpr Config(int tile_size, ClipSpaceStyle clip_space_style)
            : clip_space_style(clip_space_style), tile_size(tile_size), perspective_correction(true), use_reciprocal_for_perspective_correction(false){}
    };


    /**
 * @brief Transforms and rasterizes a CLOCKWISE triangle defined in clip space coordinates.
 *
 * This function performs per-vertex perspective division, viewport transformation,
 * and screen-space rasterization of a single triangle. The input vertices are expected
 * to be in clip space (x, y, z, w), as produced by a vertex shader. After transformation,
 * the triangle is passed to the internal rasterization stage (e.g. scanline or bounding box
 * method) defined by the configuration type.
 *
 * The rasterizer operates in tiles and invokes a user-provided fragment callback for each
 * covered region. Optional tile locking functions can be used for safe multithreaded rendering.
 *
 * ---
 *
 * @tparam CONFIG
 *     Compile-time configuration object describing rasterization behavior.
 *     It typically specifies:
 *       - `Warp` � SIMD warp type used for vectorized processing,
 *       - `clip_space_style` � clip-space convention (e.g., OpenGL, D3D),
 *       - `tile_size` � tile dimension in pixels,
 *       - `use_reciprocal_for_perspective_correction` � if true, uses a reciprocal (1/w)
 *         instruction instead of division when performing perspective-correct interpolation,
 * trading a small amount of accuracy for better performance on SIMD hardware.
 *
 * @tparam NUM_ATTRS
 *     Number of interpolated vertex attributes per vertex (e.g. UVs, color, depth, etc.).
 *
 * ---
 *
 * @param viewport
 *     Viewport description defining screen-space transformation.
 *
 * @param ax, ay, az, aw
 *     Clip-space coordinates of vertex A.
 *
 * @param bx, by, bz, bw
 *     Clip-space coordinates of vertex B.
 *
 * @param cx, cy, cz, cw
 *     Clip-space coordinates of vertex C.
 *
 * @param a_attrs
 *     Array of `NUM_ATTRS` vertex attributes for vertex A.
 *
 * @param b_attrs
 *     Array of `NUM_ATTRS` vertex attributes for vertex B.
 *
 * @param c_attrs
 *     Array of `NUM_ATTRS` vertex attributes for vertex C.
 *
 * @param fragment
 *     Fragment callback invoked for each covered region during rasterization.
 *     Signature:
 *     @code
 *     void fragment(int tx, int ty,
 *                   int column, int row,
 *                   typename CONFIG::Warp::IntMask mask,
 *                   typename CONFIG::Warp::Float* attrs);
 *     @endcode
 *     - `tx`, `ty` � tile coordinates (in tile units);
 *     - `column`, `row` � local pixel position within the tile;
 *     - `mask` � SIMD mask indicating active pixels;
 *     - `attrs` � interpolated vertex attributes for active pixels.
 *
 *     The user is responsible for writing shaded results into the framebuffer
 *     (e.g., color or depth buffers) within this callback.
 *
 * @param tile_lock
 *     Function called before processing a tile.
 *     Used for optional synchronization in multithreaded rendering.
 *
 * @param tile_unlock
 *     Function called after processing a tile.
 *     Complements `tile_lock`, used to finalize or unlock shared data.
 *
 * ---
 *
 * @note
 * - Vertex coordinates are expected to be in **clip space** before perspective division.
 * - Performs viewport transformation to convert to screen space.
 * - Supports both perspective-correct and affine interpolation depending on `CONFIG`.
 * - The rasterizer uses the top-left fill rule for consistent edge inclusion.
 * - Does not perform any implicit framebuffer writes � user logic in `fragment` controls output.
 * - If `CONFIG::perspective_correction` is enabled, the fragment shader receives an
 *       additional implicit attribute at the end of the `attrs` array: `1/w` of the fragment
 *
 * ---
 *
 * @example
 * Example usage:
 * @code
 * static constexpr Config<MyWarp> cfg {
 *     .clip_space_style = ClipSpaceStyle::OpenGL,
 *     .tile_size = 8,
 *     .perspective_correction = true
 * };
 *
 * clip_space_triangle<cfg, 3>(
 *     viewport,
 *     ax, ay, az, aw,
 *     bx, by, bz, bw,
 *     cx, cy, cz, cw,
 *     a_attrs, b_attrs, c_attrs,
 *     [&](int tx, int ty, int column, int row,
 *         MyWarp::IntMask mask, MyWarp::Float* attrs) {
 *         // Perform shading and write to color/depth buffers
 *     },
 *     [] (int tx, int ty) {}, // lock
 *     [] (int tx, int ty) {}  // unlock
 * );
 * @endcode
 */
    template<auto CONFIG, int NUM_ATTRS>
    inline constexpr void clip_space_triangle(
        Viewport viewport,
        float ax, float ay, float az, float aw,
        float bx, float by, float bz, float bw,
        float cx, float cy, float cz, float cw,
        float (&a_attrs)[NUM_ATTRS],
        float (&b_attrs)[NUM_ATTRS],
        float (&c_attrs)[NUM_ATTRS],
        auto fragment,
        auto tile_lock,
        auto tile_unlock) {

        struct ClipVertex {
            float x, y, z, w, weight0, weight1;

            inline ClipVertex lerp(const ClipVertex& other, float t) {
                return {
                    x + (other.x - x) * t,
                    y + (other.y - y) * t,
                    z + (other.z - z) * t,
                    w + (other.w - w) * t,
                    weight0 + (other.weight0 - weight0) * t,
                    weight1 + (other.weight1 - weight1) * t
                };
            }
        };

        ClipVertex vertices[10];
        int num_vertices = 3;

        vertices[0].x = ax;
        vertices[0].y = ay;
        vertices[0].z = az;
        vertices[0].w = aw;
        vertices[0].weight0 = 1;
        vertices[0].weight1 = 0;

        vertices[1].x = bx;
        vertices[1].y = by;
        vertices[1].z = bz;
        vertices[1].w = bw;
        vertices[1].weight0 = 0;
        vertices[1].weight1 = 1;

        vertices[2].x = cx;
        vertices[2].y = cy;
        vertices[2].z = cz;
        vertices[2].w = cw;
        vertices[2].weight0 = 0;
        vertices[2].weight1 = 0;

        ClipSpaceFrustum<CONFIG.clip_space_style> frustum {};

        bool all_inside = frustum.inside(vertices, num_vertices);

        if (!all_inside) {
            auto all_outside = frustum.outside(vertices, num_vertices);
            if (all_outside)
                return;

            frustum.clip_polygon(vertices, num_vertices);

            if (num_vertices < 3)
                return;
        }

        for (int i = 0; i < num_vertices; ++i) {
            auto& v = vertices[i];
            viewport.clip_space_to_screen_inv_w(v.x, v.y, v.w);
        }

        using Warp = typename decltype(CONFIG)::Warp;

        for (int i = 1; i < num_vertices - 1; ++i) {
            auto a = vertices[0];
            auto b = vertices[i];
            auto c = vertices[i + 1];

             if constexpr (CONFIG.perspective_correction) {

                 float new_a_attrs[NUM_ATTRS + 1];
                 float new_b_attrs[NUM_ATTRS + 1];
                 float new_c_attrs[NUM_ATTRS + 1];

                 for (int j = 0; j < NUM_ATTRS; ++j) {
                     new_a_attrs[j] = (a_attrs[j] * a.weight0 + b_attrs[j] * a.weight1 + c_attrs[j] * (1.f - a.weight0 - a.weight1)) * a.w;
                     new_b_attrs[j] = (a_attrs[j] * b.weight0 + b_attrs[j] * b.weight1 + c_attrs[j] * (1.f - b.weight0 - b.weight1)) * b.w;
                     new_c_attrs[j] = (a_attrs[j] * c.weight0 + b_attrs[j] * c.weight1 + c_attrs[j] * (1.f - c.weight0 - c.weight1)) * c.w;
                 }

                 constexpr int INV_W = NUM_ATTRS;

                 new_a_attrs[INV_W] = a.w;
                 new_b_attrs[INV_W] = b.w;
                 new_c_attrs[INV_W] = c.w;

                 auto perspective_correct_fragment = [&](int tx, int ty, int column, int row, typename Warp::IntMask mask, typename Warp::Float* attrs) {

                     typename Warp::Float attrs_copy[NUM_ATTRS + 1];

                     for (int i = 0; i < NUM_ATTRS; ++i)
                         if constexpr (CONFIG.use_reciprocal_for_perspective_correction)
                             attrs_copy[i] = Warp::mul(attrs[i], Warp::reciprocal(attrs[INV_W]));
                         else
                             attrs_copy[i] = Warp::div(attrs[i], attrs[INV_W]);

                     attrs_copy[INV_W] = attrs[INV_W];

                     fragment(tx, ty, column, row, mask, attrs_copy);
                 };

                 triangle<Warp, CONFIG.tile_size>(
                     a.x, a.y,
                     b.x, b.y,
                     c.x, c.y,
                     new_a_attrs,
                     new_b_attrs,
                     new_c_attrs,
                     perspective_correct_fragment,
                     tile_lock,
                     tile_unlock);
             }
             else 
            {
                float new_a_attrs[NUM_ATTRS];
                float new_b_attrs[NUM_ATTRS];
                float new_c_attrs[NUM_ATTRS];

                for (int j = 0; j < NUM_ATTRS; ++j) {
                    new_a_attrs[j] = (a_attrs[j] * a.weight0 + b_attrs[j] * a.weight1 + c_attrs[j] * (1.f - a.weight0 - a.weight1));
                    new_b_attrs[j] = (a_attrs[j] * b.weight0 + b_attrs[j] * b.weight1 + c_attrs[j] * (1.f - b.weight0 - b.weight1));
                    new_c_attrs[j] = (a_attrs[j] * c.weight0 + b_attrs[j] * c.weight1 + c_attrs[j] * (1.f - c.weight0 - c.weight1));
                }
                
                triangle<Warp, CONFIG.tile_size>(
                    a.x, a.y,
                    b.x, b.y,
                    c.x, c.y,
                    new_a_attrs,
                    new_b_attrs,
                    new_c_attrs,
                    fragment,
                    tile_lock,
                    tile_unlock);
            }
        }
    }
}