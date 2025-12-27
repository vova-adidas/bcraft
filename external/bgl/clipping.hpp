#pragma once

namespace bgl {
  
    enum class ClipSpaceStyle {
        DirectX,
        OpenGL
    };

    enum class ClipPlane {
        Near,
        Far,
        Left,
        Right,
        Bot,
        Top
    };

    template<ClipSpaceStyle Style>
    class ClipSpaceFrustum {
    private:

        template <typename A, typename B>
        struct Pair {
            A inside_test;
            B solution;
        };

        template<typename A, typename B>
        Pair(A, B) -> Pair<A, B>;

        template<ClipPlane PLANE>
        constexpr static auto plane_funcs() {

            auto near_solution = [](auto a, auto b) {
                if constexpr (Style == ClipSpaceStyle::DirectX) {
                    auto res = a.lerp(b, (0 - a.z) / (b.z - a.z));
                    //to fix float interpolation error
                    res.z = 0;
                    return res;
                }
                else {
                    auto res = a.lerp(b, (-a.w - a.z) / ((b.z - a.z) + (b.w - a.w)));
                    //to fix float interpolation error
                    res.z = -a.w;
                    return res;
                }
                };

            auto far_solution = [](auto a, auto b) {
                auto res = a.lerp(b, (a.w - a.z) / ((b.z - a.z) - (b.w - a.w)));
                //to fix float interpolation error
                res.z = (res.z >= 0) ? res.z : -res.z;
                return res;
                };
            auto left_solution = [](auto a, auto b) {
                auto res = a.lerp(b, (-a.w - a.x) / ((b.x - a.x) + (b.w - a.w)));
                //to fix float interpolation error
                res.x = (res.x >= 0) ? res.w : -res.w;
                return res;
                };
            auto right_solution = [](auto a, auto b) {
                auto res = a.lerp(b, (a.w - a.x) / ((b.x - a.x) - (b.w - a.w)));
                //to fix float interpolation error
                res.x = (res.x >= 0) ? res.w : -res.w;
                return res;
                };
            auto bot_solution = [](auto a, auto b) {
                auto res = a.lerp(b, (-a.w - a.y) / ((b.y - a.y) + (b.w - a.w)));
                //to fix float interpolation error
                res.y = (res.y >= 0) ? res.w : -res.w;
                return res;
                };
            auto top_solution = [](auto a, auto b) {
                auto res = a.lerp(b, (a.w - a.y) / ((b.y - a.y) - (b.w - a.w)));
                //to fix float interpolation error
                res.y = (res.y >= 0) ? res.w : -res.w;
                return res;
                };

            if constexpr (PLANE == ClipPlane::Near) {
                return Pair { 
                    [](auto px, auto py, auto pz, auto pw) { 
                        if constexpr (Style == ClipSpaceStyle::DirectX) {
                            return pz >= 0;
                        }
                        else
                            return pz >= -pw; 
                        }, near_solution };
            }
            else if constexpr (PLANE == ClipPlane::Far)
                return Pair{ [](auto px, auto py, auto pz, auto pw) { return pz <= pw; }, far_solution };
            else if constexpr (PLANE == ClipPlane::Left)
                return Pair{ [](auto px, auto py, auto pz, auto pw) { return px >= -pw; }, left_solution };
            else if constexpr (PLANE == ClipPlane::Right)
                return Pair{ [](auto px, auto py, auto pz, auto pw) { return px <= pw; }, right_solution };
            else if constexpr (PLANE == ClipPlane::Bot)
                return Pair{ [](auto px, auto py, auto pz, auto pw) { return py >= -pw; }, bot_solution };
            else
                return Pair{ [](auto px, auto py, auto pz, auto pw) { return py <= pw; }, top_solution };
        }
    public:

        template<ClipPlane PLANE>
        inline constexpr static bool inside_plane(float px, float py, float pz, float pw) { return plane_funcs<PLANE>().inside_test(px, py, pz, pw); };

        inline constexpr static bool inside(float px, float py, float pz, float pw) {

            return inside_plane<ClipPlane::Near>(px, py, pz, pw)
                &&
                inside_plane<ClipPlane::Far>(px, py, pz, pw)
                &&
                inside_plane<ClipPlane::Left>(px, py, pz, pw)
                &&
                inside_plane<ClipPlane::Right>(px, py, pz, pw)
                &&
                inside_plane<ClipPlane::Bot>(px, py, pz, pw)
                &&
                inside_plane<ClipPlane::Top>(px, py, pz, pw);
        }

        /**
         * @brief Checks whether all vertices of a polygon are inside the view frustum.
         *
         * @tparam TVertex
         *     Vertex type. Must provide public members `x`, `y`, `z`, `w` representing
         *     homogeneous clip-space coordinates.
         *
         * @param vertices
         *     Pointer to an array of `TVertex` representing the polygon vertices.
         *
         * @param num_vertices
         *     Number of vertices in the polygon.
         *
         * @return
         *     `true` if all vertices are inside the clip space; `false` otherwise.
         */
        template<typename TVertex>
        inline constexpr static bool inside(TVertex* vertices, int num_vertices) {
            bool all_inside = true;
            for (int i = 0; i < num_vertices; ++i) {
                if (!inside(vertices[i].x, vertices[i].y, vertices[i].z, vertices[i].w)) {
                    all_inside = false;
                    break;
                }
            }

            return all_inside;
        }

        /**
         * @brief Checks whether all vertices of a polygon are completely outside the view frustum.
         *
         * @tparam TVertex
         *     Vertex type. Must provide public members `x`, `y`, `z`, `w` representing
         *     homogeneous clip-space coordinates.
         *
         * @param vertices
         *     Pointer to an array of `TVertex` representing the polygon vertices.
         *
         * @param num_vertices
         *     Number of vertices in the polygon.
         *
         * @return
         *     `true` if all vertices are inside the clip space; `false` otherwise.
         */
        template<typename TVertex>
        inline constexpr static bool outside(TVertex* vertices, int num_vertices) {
            bool all_outside = true;
            for (int i = 0; i < num_vertices; ++i) {
                if (inside_plane<ClipPlane::Near>(vertices[i].x, vertices[i].y, vertices[i].z, vertices[i].w)) {
                    all_outside = false;
                    break;
                }
            }
            if (all_outside) return true;

            all_outside = true;
            for (int i = 0; i < num_vertices; ++i) {
                if (inside_plane<ClipPlane::Far>(vertices[i].x, vertices[i].y, vertices[i].z, vertices[i].w)) {
                    all_outside = false;
                    break;
                }
            }
            if (all_outside) return true;

            all_outside = true;
            for (int i = 0; i < num_vertices; ++i) {
                if (inside_plane<ClipPlane::Left>(vertices[i].x, vertices[i].y, vertices[i].z, vertices[i].w)) {
                    all_outside = false;
                    break;
                }
            }
            if (all_outside) return true;

            all_outside = true;
            for (int i = 0; i < num_vertices; ++i) {
                if (inside_plane<ClipPlane::Right>(vertices[i].x, vertices[i].y, vertices[i].z, vertices[i].w)) {
                    all_outside = false;
                    break;
                }
            }
            if (all_outside) return true;

            all_outside = true;
            for (int i = 0; i < num_vertices; ++i) {
                if (inside_plane<ClipPlane::Bot>(vertices[i].x, vertices[i].y, vertices[i].z, vertices[i].w)) {
                    all_outside = false;
                    break;
                }
            }
            if (all_outside) return true;

            all_outside = true;
            for (int i = 0; i < num_vertices; ++i) {
                if (inside_plane<ClipPlane::Top>(vertices[i].x, vertices[i].y, vertices[i].z, vertices[i].w)) {
                    all_outside = false;
                    break;
                }
            }
            if (all_outside) return true;

            return false;
        }


        /**
         * @brief Clips a convex polygon in homogeneous clip space in-place.
         *
         * This function clips a convex polygon against the view plane
         * directly in the input array. The number of vertices is updated to reflect the clipped polygon.
         * Additional vertex attributes can be interpolated using `TVertex.lerp`.
         *
         * ---
         *
         * @tparam TVertex
         *     Vertex type. Must provide:
         *       - Public floating-point members `x`, `y`, `z`, `w` (homogeneous clip-space coordinates).
         *       - A method:
         *         @code
         *         TVertex lerp(const TVertex& v, float t);
         *         @endcode
         *         used to interpolate all additional attributes when generating new vertices.
         *
         * @param vertices
         *     Pointer to an array of `TVertex` representing the polygon vertices.
         *     The array is modified in-place to store the clipped polygon.
         *     **Must be large enough to hold all vertices generated during clipping.**
         *
         * @param num_vertices
         *     Input: number of vertices in the polygon.
         *     Output: number of vertices after clipping.
         *
         * ---
         */
        template<ClipPlane PLANE, typename TVertex>
        inline constexpr static void clip_polygon_by_plane(TVertex* vertices, int& num_vertices) {
 
            int v_index = 0;

            const auto add_vertex = [&](auto v) { vertices[v_index] = v; v_index++;};
               
            auto pf = plane_funcs<PLANE>();

            const auto clip_edge = [&](auto a, auto b) {
                auto a_inside = pf.inside_test(a.x, a.y, a.z, a.w);
                auto b_inside = pf.inside_test(b.x, b.y, b.z, b.w);
                if (a_inside ^ b_inside) {
                    add_vertex(pf.solution(a, b));

                }
                if (b_inside)
                    add_vertex(b);

                };

            auto a = vertices[0];
            auto first = a;

            TVertex b{};
            for (int i = 1; i < num_vertices; ++i) {
                b = vertices[i];
                clip_edge(a, b);
                a = b;
            }

            b = first;
            clip_edge(a, b);

            num_vertices = v_index;
        }

       /**
        * @brief Clips a convex polygon in homogeneous clip space in-place.
        *
        * This function clips a convex  polygon against the frustum
        * directly in the input array. The number of vertices is updated to reflect the clipped polygon.
        * Additional vertex attributes can be interpolated using `TVertex.lerp`.
        *
        * ---
        *
        * @tparam TVertex
        *     Vertex type. Must provide:
        *       - Public floating-point members `x`, `y`, `z`, `w` (homogeneous clip-space coordinates).
        *       - A method:
        *         @code
        *         TVertex lerp(const TVertex& v, float t);
        *         @endcode
        *         used to interpolate all additional attributes when generating new vertices.
        *
        * @param vertices
        *     Pointer to an array of `TVertex` representing the polygon vertices.
        *     The array is modified in-place to store the clipped polygon.
        *     **Must be large enough to hold all vertices generated during clipping.**
        *
        * @param num_vertices
        *     Input: number of vertices in the polygon.
        *     Output: number of vertices after clipping.
        *
        * ---
        * ---
        *
        * @example
        * @code
        * struct Vertex {
        *     float x, y, z, w;
        *     float u, v; // extra attributes
        *
        *     static Vertex lerp(const Vertex& o, float t) {
        *         Vertex v;
        *         o.x = x + t*(o.x - x);
        *         ...
        *         v.v = v + t*(o.v - v);
        *         return v;
        *     }
        * };
        *
        * Vertex poly[16]; // make sure this array is large enough for possible new vertices
        * int n = 6;
        * frustum.clip_polygon(poly, n);
        * // poly now contains n vertices after clipping
        * @endcode
        */
        template<typename TVertex>
        inline constexpr static void clip_polygon(TVertex* vertices, int& num_vertices) {

            auto need_clip_near = false;
            auto need_clip_far = false;
            auto need_clip_left = false;
            auto need_clip_right = false;
            auto need_clip_bot = false;
            auto need_clip_top = false;

            for (int i = 0; i < num_vertices; ++i) {
                if (!inside_plane<ClipPlane::Near>(vertices[i].x, vertices[i].y, vertices[i].z, vertices[i].w))
                    need_clip_near = true;

                if (!inside_plane<ClipPlane::Far>(vertices[i].x, vertices[i].y, vertices[i].z, vertices[i].w))
                    need_clip_far = true;

                if (!inside_plane<ClipPlane::Left>(vertices[i].x, vertices[i].y, vertices[i].z, vertices[i].w))
                    need_clip_left = true;

                if (!inside_plane<ClipPlane::Right>(vertices[i].x, vertices[i].y, vertices[i].z, vertices[i].w))
                    need_clip_right = true;

                if (!inside_plane<ClipPlane::Bot>(vertices[i].x, vertices[i].y, vertices[i].z, vertices[i].w))
                    need_clip_bot = true;

                if (!inside_plane<ClipPlane::Top>(vertices[i].x, vertices[i].y, vertices[i].z, vertices[i].w))
                    need_clip_top = true;
            }

            if (need_clip_near)
                clip_polygon_by_plane<ClipPlane::Near>(vertices, num_vertices);

            if (need_clip_far)
                clip_polygon_by_plane<ClipPlane::Far>(vertices, num_vertices);

            if (need_clip_left)
                clip_polygon_by_plane<ClipPlane::Left>(vertices, num_vertices);

            if (need_clip_right)
                clip_polygon_by_plane<ClipPlane::Right>(vertices, num_vertices);

            if (need_clip_bot)
                clip_polygon_by_plane<ClipPlane::Bot>(vertices, num_vertices);

            if (need_clip_top)
                clip_polygon_by_plane<ClipPlane::Top>(vertices, num_vertices);
        }
    };
}