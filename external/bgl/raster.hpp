#pragma once
#include <limits>

namespace bgl {
    
    /**
    * @class TriangleScanlineIterator
    * @brief Iterates over all scanlines of a triangle, covering every pixel
    *        that is touched by the triangle edges, including partial coverage.
    *
    * This class generates a series of scanlines (horizontal lines) that
    * intersect the triangle. For each scanline, it provides the range
    * of X coordinates that the triangle covers. It guarantees that all
    * pixels that the triangle touches are included.
    *
    * Features:
    *  - Full coverage of edge pixels.
    *  - Handles triangles of arbitrary orientation.
    *
    * Usage:
    * @code
    * TriangleScanlineIterator iter(ax, ay, bx, by, cx, cy);
    * while (iter.next()) {
    *     for (int x = iter.begin(); x < iter.end(); ++x) {
    *         int y = iter.line();
    *         // process pixel (x, y)
    *     }
    * }
    * @endcode
    */

    template<typename Float>
    class TriangleScanlineIteratorBase {
    public:
        constexpr static Float Epsilon = 0.000001;

        constexpr static auto swap(auto& a, auto& b) {
            auto temp = a;
            a = b;
            b = temp;
        };

        constexpr static auto ceil(Float x) {
            int i = static_cast<int>(x);
            return x > i ? i + 1 : i;
        };

        constexpr static auto floor(Float x) {
            int i = static_cast<int>(x);
            return i;
        };

        constexpr static auto float_eq(auto x, auto y) {
            auto d = x - y;
            if (d < 0)
                d = -d;
            return d < Epsilon;
        };

        Float
            left[6],
            right[6],
            delta_left[5]{},
            delta_right[5]{},
            cur_x0,
            cur_x1, 
            cur_begin,
            cur_end;

        int 
            n[5], 
            k = 0, 
            i = 0,
            j = 0,
            cur_y,
            cur_line;

    public:
        constexpr TriangleScanlineIteratorBase(const TriangleScanlineIteratorBase&) = delete;

        constexpr TriangleScanlineIteratorBase& operator=(const TriangleScanlineIteratorBase&) = delete;

        constexpr TriangleScanlineIteratorBase(Float ax, Float ay, Float bx, Float by, Float cx, Float cy) {
            
            if (by > cy) {
                swap(bx, cx);
                swap(by, cy);
            }

            if (ay > by) {
                swap(ax, bx);
                swap(ay, by);
            }

            if (by > cy) {
                swap(bx, cx);
                swap(by, cy);
            }

            const auto saturate = [](Float x) { return x < 0 ? 0 : (x < 1 ? x : 1); };

            int min_y = floor(ay);
            int center_y = floor(by);
            int max_y = ceil(cy);
            
            if (static_cast<int>(ay) == static_cast<int>(cy)) {
                auto x0 = (ax < bx ? (ax < cx ? ax : cx) : (bx < cx ? bx : cx));
                auto x1 = (ax > bx ? (ax > cx ? ax : cx) : (bx > cx ? bx : cx));

                left[0] = x0;
                right[0] = x1;
                n[0] = 1;
                for (int i = 1; i < 5; ++i)
                    n[i] = 0;
                k++;
            }
            else
            {

                auto mid_y = by;
                auto left_x = bx;
                auto right_x = ax + (cx - ax) * (mid_y - ay) / (cy - ay);

                if (left_x > right_x)
                    swap(left_x, right_x);

                auto left_edge_correction = left_x - ax < 0 ? 1 : 0;
                auto right_edge_correction = right_x - ax > 0 ? 1 : 0;

                
                auto inv_by_ay = 1 / (by - ay);
                if (min_y < center_y) {
                    auto y = min_y;
                    auto t0 = (y + left_edge_correction - ay) * inv_by_ay;
                    auto t1 = (y + right_edge_correction - ay) * inv_by_ay;

                    left[k] = ax + (left_x - ax) * saturate(t0);
                    right[k] = (ax + (right_x - ax) * saturate(t1));
                    n[k] = 1;
                    k++;
                }

                if (center_y - (min_y + 1) > 0) {
                    
                    n[k] = center_y - (min_y + 1);

                    {
                        auto init = ax + (right_x - ax) * ((min_y + 1 - ay) * inv_by_ay);
                        auto end = ax + (right_x - ax) * ((center_y - ay) * inv_by_ay);

                        right[k] = init;
                        delta_right[k] = (end - init) / n[k];
                        right[k] += right_edge_correction * delta_right[k];
                    }

                    {
                        auto init = ax + (left_x - ax) * ((min_y + 1 - ay) * inv_by_ay);
                        auto end = ax + (left_x - ax) * ((center_y - ay) * inv_by_ay);

                        left[k] = init;
                        delta_left[k] = (end - init) / n[k];
                        left[k] += left_edge_correction * delta_left[k];
                    }
                    k++;
                }

                auto inv_cy_mid_y = 1 / (cy - mid_y);

                {
                    auto y = center_y;

                    Float x0 = left_x;
                    auto edge = float_eq(by, ay) ? left_x : (ax + (left_x - ax) * saturate((y - ay) * inv_by_ay));
                    x0 = edge < x0 ? edge : x0;
                    edge = float_eq(cy, mid_y) ? left_x : (left_x + (cx - left_x) * saturate((y + 1 - mid_y) * inv_cy_mid_y));
                    x0 = edge < x0 ? edge : x0;

                    Float x1 = right_x;
                    edge = float_eq(by, ay) ? right_x : (ax + (right_x - ax) * saturate((y - ay) * inv_by_ay));
                    x1 = edge > x1 ? edge : x1;
                    edge = float_eq(cy, mid_y) ? right_x : (right_x + (cx - right_x) * saturate((y + 1 - mid_y) * inv_cy_mid_y));
                    x1 = edge > x1 ? edge : x1;

                    n[k] = 1;
                    left[k] = x0;
                    right[k] = x1;
                    ++k;
                }

                left_edge_correction = cx - left_x < 0 ? 1 : 0;
                right_edge_correction = cx - right_x > 0 ? 1 : 0;

                if(max_y - 1 - (center_y + 1) > 0) {
                    n[k] = max_y - 1 - (center_y + 1);
                    {
                        auto init = right_x + (cx - right_x) * ((center_y + 1 - mid_y) * inv_cy_mid_y);
                        auto end = right_x + (cx - right_x) * ((max_y - 1 - mid_y) * inv_cy_mid_y);

                        right[k] = init;
                        delta_right[k] = (end - init) / n[k];
                        right[k] += right_edge_correction * delta_right[k];
                    }

                    {
                        auto init = left_x + (cx - left_x) * ((center_y + 1 - mid_y) * inv_cy_mid_y);
                        auto end = left_x + (cx - left_x) * ((max_y - 1 - mid_y) * inv_cy_mid_y);

                        left[k] = init;
                        delta_left[k] = (end - init) / n[k];
                        left[k] += left_edge_correction * delta_left[k];
                    }
                    ++k;
                }

                {
                    auto y = max_y - 1;
                    if (center_y < y) {
                        auto t0 = (y + left_edge_correction - mid_y) * inv_cy_mid_y;
                        auto t1 = (y + right_edge_correction - mid_y) * inv_cy_mid_y;

                        n[k] = 1;
                        left[k] = ((left_x + (cx - left_x) * saturate(t0)));
                        right[k] = ((right_x + (cx - right_x) * saturate(t1)));
                        ++k;
                    }
                }
            }

            cur_y = min_y;
            cur_x0 = left[j];
            cur_x1 = right[j];
            cur_line = cur_y;
        }

        constexpr bool next() {
            
            if (j >= k)
                return false;
       
            cur_begin = cur_x0;
            cur_end = cur_x1;
            cur_line = cur_y;
            
            cur_x0 += delta_left[j];
            cur_x1 += delta_right[j];

            ++cur_y;
            ++i;


            if (i == n[j]) {
                ++j;
                i = 0;

                cur_x0 = left[j];
                cur_x1 = right[j];
            }
           

            return true;
        }

        constexpr Float begin() const { return cur_begin; }

        constexpr Float end() const { return cur_end; }

        constexpr int line() const { return cur_line; }
    };

    using TriangleScanlineIterator = TriangleScanlineIteratorBase<float>;

    template<typename Float>
    class EllipseScanlineIteratorBase {
    public:
        constexpr static Float Epsilon = 0.000001;

        // Ваша реализация floor
        constexpr static int floor(Float x) {
            int i = static_cast<int>(x);
            return (x < i) ? i - 1 : i;
        }

        // Ваша реализация ceil
        constexpr static int ceil(Float x) {
            int i = static_cast<int>(x);
            return (x > i) ? i + 1 : i;
        }

        // Собственная реализация SQRT (метод Ньютона) для constexpr или runtime без math.h
        constexpr static Float sqrt(Float x) {
            if (x <= 0) return 0;
            Float res = x;
            // 5-7 итераций достаточно для точности растеризации
            for (int i = 0; i < 6; ++i) {
                res = 0.5f * (res + x / res);
            }
            return res;
        }

    private:
        Float cx, cy, rx, ry;
        Float inv_ry_sq; // Предрассчитанное значение 1 / (ry * ry)
        int min_y, max_y;
        int cur_y;

        Float cur_begin;
        Float cur_end;

    public:
        constexpr EllipseScanlineIteratorBase(const EllipseScanlineIteratorBase&) = delete;

        constexpr EllipseScanlineIteratorBase(Float centerX, Float centerY, Float radiusX, Float radiusY)
            : cx(centerX), cy(centerY), rx(radiusX), ry(radiusY)
        {
            min_y = floor(cy - ry);
            max_y = ceil(cy + ry);
            cur_y = min_y;

            // Предрассчитываем инверсию квадрата радиуса для ускорения цикла
            inv_ry_sq = 1.0f / (ry * ry);
        }

        constexpr bool next() {
            if (cur_y >= max_y)
                return false;

            Float dy = (static_cast<Float>(cur_y) + 0.5f) - cy;
            Float dy_sq = dy * dy;


            Float term = 1.0f - (dy_sq * inv_ry_sq);

            if (term < 0) {
                cur_begin = cx;
                cur_end = cx;
            } else {
                Float width = rx * sqrt(term);
                cur_begin = cx - width;
                cur_end = cx + width;
            }

            cur_y++;
            return true;
        }

        constexpr Float begin() const { return cur_begin; }
        constexpr Float end() const { return cur_end; }
        constexpr int line() const { return cur_y - 1; }
    };
    using EllipseScanlineIterator = EllipseScanlineIteratorBase<float>;

    //Edgefunc > 0 or >= 0
    template<typename TWarp, int SubpixelRes>
    struct EdgeInequalityIterator {
        constexpr static int LIMIT_MIN = -(1 << 28);
        constexpr static int LIMIT_MAX = (1 << 28);
        using bigint = long long;

        int a{}, q{}, r{}, dq{}, dr{};
        bool flip{}, horizontal{};

        typename TWarp::Int vec_q, vec_r;

        constexpr EdgeInequalityIterator() = default;

        constexpr EdgeInequalityIterator(bool strict, bigint  ax, bigint ay, bigint  bx, bigint  by, bigint init_y) {


            a = by - ay;

            bigint c = ax - bx;

            bigint non_strict_correction = strict ? 0 : 1;

            long long b64 = c * init_y + static_cast<bigint>(bx) * ay - static_cast<bigint>(ax) * by + non_strict_correction;

            flip = (a < 0);
            this->a = flip ? -a : a;

            bigint dq1, dr1;
            bigint q64, r64, dq64, dr64, dq164, dr164, _;

            const auto euclid_div = [](bigint a, bigint b, bigint& q, bigint& r) {

                q = a / b;
                r = a % b;

                if (r < 0) {
                    q -= 1;
                    r += b;
                }
            };

            const auto horizontal_setup = [&]() {
                horizontal = true;
                a = LIMIT_MAX;

                if (c > 0)
                {
                    q64 = init_y - ay;
                    dq64 = SubpixelRes * TWarp::SIZE;
                    dq164 = SubpixelRes;
                }
                else
                {
                    q64 = ay - init_y;
                    dq64 = -SubpixelRes * TWarp::SIZE;
                    dq164 = -SubpixelRes;

                }

                r64 = 0;
                dr64 = 0;
                dr164 = 0;
            };

            const auto default_setup = [&]() {
                horizontal = false;

                euclid_div(-b64, SubpixelRes, q64, _);
                euclid_div(q64, this->a, q64, r64);

                if (q64 >= LIMIT_MAX || q64 <= LIMIT_MIN)
                    horizontal_setup();
                else {
                    euclid_div(-c * TWarp::SIZE, this->a, dq64, dr64);
                    euclid_div(-c, this->a, dq164, dr164);
                    q64++;
                }
            };

            if (this->a != 0)
                default_setup();
            else
                horizontal_setup();

            q = static_cast<int>(q64);
            r = static_cast<int>(r64);
            dq = static_cast<int>(dq64);
            dr = static_cast<int>(dr64);
            dq1 = static_cast<int>(dq164);
            dr1 = static_cast<int>(dr164);

            alignas(TWarp::ALIGNMENT) int a_vec_q[TWarp::SIZE];
            alignas(TWarp::ALIGNMENT) int a_vec_r[TWarp::SIZE];
            int r1 = r;
            int q1 = q;

            for (int i = 0; i < TWarp::SIZE; ++i) {
                a_vec_q[i] = q1;
                a_vec_r[i] = r1;

                r1 += dr1;
                if (r1 >= a) {
                    ++q1;
                    r1 -= a;
                }
                q1 += dq1;
            }
            vec_q = TWarp::load(a_vec_q);
            vec_r = TWarp::load(a_vec_r);
        }

        constexpr void next_intervals(typename TWarp::Int& min, typename TWarp::Int& max) {
            auto vec_dr = TWarp::scalar(dr);
            auto vec_dq = TWarp::scalar(dq);
            auto vec_a = TWarp::scalar(a);
            auto neg_vec_q = TWarp::neg(vec_q);
            auto one = TWarp::scalar(1);
            auto zero = TWarp::scalar(0);

            auto vec_limit_min = TWarp::scalar(LIMIT_MIN);
            auto vec_limit_max = TWarp::scalar(LIMIT_MAX);

            auto h_x_min = TWarp::select(TWarp::greater(vec_q, zero), vec_limit_min, vec_limit_max);
            auto h_x_max = TWarp::select(TWarp::greater(vec_q, zero), vec_limit_max, vec_limit_min);

            auto x_min = TWarp::select(TWarp::mask_from_bool(flip), vec_limit_min, vec_q);
            auto x_max = TWarp::select(TWarp::mask_from_bool(flip), neg_vec_q, vec_limit_max);

            x_min = horizontal ? h_x_min : x_min;
            x_max = horizontal ? h_x_max : x_max;

            vec_r = TWarp::add(vec_r, vec_dr);

            vec_q = TWarp::select(TWarp::greater_equals(vec_r, vec_a), vec_q + one, vec_q);
            vec_r = TWarp::select(TWarp::greater_equals(vec_r, vec_a), vec_r - vec_a, vec_r);

            vec_q = TWarp::add(vec_q, vec_dq);

            min = x_min;
            max = x_max;
        }
    };

    template<typename TWarp, int TILE_SIZE, int NUM_ATTRS>
    void triangle(
        float fax,
        float fay,
        float fbx,
        float fby,
        float fcx,
        float fcy,
        float(&a_attributes)[NUM_ATTRS],
        float(&b_attributes)[NUM_ATTRS],
        float(&c_attributes)[NUM_ATTRS],
        auto fragment,
        auto tile_lock,
        auto tile_unlock) {

        const auto div_floor = [](long long a, long long b) {
            long long q = a / b;
            long long r = a % b;
            if (r != 0 && ((r > 0) != (b > 0)))
                q -= 1;
            return q;
        };

        const auto div_ceil = [](long long a, long long b) {
            long long q = a / b;
            long long r = a % b;
            if (r != 0 && ((r > 0) == (b > 0)))
                q += 1;
            return q;
        };

        constexpr int SubpixelRes = 1 << 11;

        alignas(TWarp::ALIGNMENT) int temp_idx[TWarp::SIZE];
        alignas(TWarp::ALIGNMENT) float temp_fidx[TWarp::SIZE];

        for (int i = 0; i < TWarp::SIZE; ++i) {
            temp_idx[i] = i;
            temp_fidx[i] = static_cast<float>(i);
        }

        auto idx = TWarp::load(temp_idx);
        auto fidx = TWarp::load(temp_fidx);

        auto ax = static_cast<long long>(static_cast<double>(fax) * SubpixelRes + .5);
        auto ay = static_cast<long long>(static_cast<double>(fay) * SubpixelRes + .5);
        auto bx = static_cast<long long>(static_cast<double>(fbx) * SubpixelRes + .5);
        auto by = static_cast<long long>(static_cast<double>(fby) * SubpixelRes + .5);
        auto cx = static_cast<long long>(static_cast<double>(fcx) * SubpixelRes + .5);
        auto cy = static_cast<long long>(static_cast<double>(fcy) * SubpixelRes + .5);

        auto edge = [](auto ax, auto ay, auto bx, auto by, auto x, auto y) {
            return (by - ay) * x + (ax - bx) * y + (bx * ay - ax * by);
        };

        long long area = edge(ax, ay, bx, by, cx, cy);
        if (area <= 0)
            return;

        double qax = (1. * ax / SubpixelRes);
        double qbx = (1. * bx / SubpixelRes);
        double qcx = (1. * cx / SubpixelRes);
        double qay = (1. * ay / SubpixelRes);
        double qby = (1. * by / SubpixelRes);
        double qcy = (1. * cy / SubpixelRes);

        long long min_y = ay < by ? (ay < cy ? ay : cy) : (by < cy ? by : cy);
        long long max_y = ay > by ? (ay > cy ? ay : cy) : (by > cy ? by : cy);

        long long ty_start = div_floor(min_y, (long long) SubpixelRes * TILE_SIZE);
        long long ty_end = div_ceil(max_y, (long long) SubpixelRes * TILE_SIZE);

        ax -= SubpixelRes / 2;
        ay -= SubpixelRes / 2;
        bx -= SubpixelRes / 2;
        by -= SubpixelRes / 2;
        cx -= SubpixelRes / 2;
        cy -= SubpixelRes / 2;

        auto is_top_left = [](long long ax, long long ay, long long bx, long long by) {
            long long dx = (bx - ax);
            long long dy = (by - ay);
            return (dy > 0) || (dy == 0 && dx < 0);
        };
        
        auto edge0 = EdgeInequalityIterator<TWarp, SubpixelRes>(!is_top_left(bx, by, cx, cy), bx, by, cx, cy, (long long)ty_start * TILE_SIZE * SubpixelRes);
        auto edge1 = EdgeInequalityIterator<TWarp, SubpixelRes>(!is_top_left(cx, cy, ax, ay), cx, cy, ax, ay, (long long)ty_start * TILE_SIZE * SubpixelRes);
        auto edge2 = EdgeInequalityIterator<TWarp, SubpixelRes>(!is_top_left(ax, ay, bx, by), ax, ay, bx, by, (long long)ty_start * TILE_SIZE * SubpixelRes);

        double scalar_start[3]{
            edge(qbx, qby, qcx, qcy, 0.5, 0.5),
            edge(qcx, qcy, qax, qay, 0.5, 0.5),
            edge(qax, qay, qbx, qby, 0.5, 0.5)
        };

        double inv_area = 1.0 / edge(qax, qay, qbx, qby, qcx, qcy);

        double scalar_dx[3] {
            (qcy - qby),
            (qay - qcy),
            (qby - qay)
        };

        double scalar_dy[3] {
            (qbx - qcx),
            (qcx - qax),
            (qax - qbx)
        };

        alignas(TWarp::ALIGNMENT) typename TWarp::Float a_start[NUM_ATTRS];
        alignas(TWarp::ALIGNMENT) typename TWarp::Float a_dx[NUM_ATTRS];
        alignas(TWarp::ALIGNMENT) typename TWarp::Float a_dy[NUM_ATTRS];
        alignas(TWarp::ALIGNMENT) typename TWarp::Float a_dx_wide[NUM_ATTRS];

        auto w0_dx = TWarp::scalar(static_cast<float>(scalar_dx[0] * inv_area));
        auto w1_dx = TWarp::scalar(static_cast<float>(scalar_dx[1] * inv_area));
        auto w2_dx = TWarp::scalar(static_cast<float>(scalar_dx[2] * inv_area));
        auto w0_dy = TWarp::scalar(static_cast<float>(scalar_dy[0] * inv_area));
        auto w1_dy = TWarp::scalar(static_cast<float>(scalar_dy[1] * inv_area));
        auto w2_dy = TWarp::scalar(static_cast<float>(scalar_dy[2] * inv_area));

        auto w0_start = TWarp::fma(fidx, w0_dx, TWarp::scalar(static_cast<float>(scalar_start[0] * inv_area)));
        auto w1_start = TWarp::fma(fidx, w1_dx, TWarp::scalar(static_cast<float>(scalar_start[1] * inv_area)));
        auto w2_start = TWarp::fma(fidx, w2_dx, TWarp::scalar(static_cast<float>(scalar_start[2] * inv_area)));

        for (int i = 0; i < NUM_ATTRS; ++i) {
            auto attr_a = TWarp::scalar(a_attributes[i]);
            auto attr_b = TWarp::scalar(b_attributes[i]);
            auto attr_c = TWarp::scalar(c_attributes[i]);

            a_start[i] = TWarp::fma(attr_a, w0_start, TWarp::fma(attr_b, w1_start, TWarp::mul(attr_c, w2_start)));
            a_dx[i] = TWarp::fma(attr_a, w0_dx, TWarp::fma(attr_b, w1_dx, TWarp::mul(attr_c, w2_dx)));
            a_dy[i] = TWarp::fma(attr_a, w0_dy, TWarp::fma(attr_b, w1_dy, TWarp::mul(attr_c, w2_dy)));
            a_dx_wide[i] = TWarp::mul(a_dx[i], TWarp::scalar(TWarp::SIZE * 1.f));
        }

        static constexpr int BATCH_SIZE = 8;

        int t_left[BATCH_SIZE];
        int t_right[BATCH_SIZE];
        alignas(TWarp::ALIGNMENT) int left[BATCH_SIZE * TILE_SIZE];
        alignas(TWarp::ALIGNMENT) int right[BATCH_SIZE * TILE_SIZE];

        for (int ty = ty_start; ty < ty_end;) {
            int b = 0;

            while (b != BATCH_SIZE) {
                constexpr int LimitMin = EdgeInequalityIterator<TWarp, SubpixelRes>::LIMIT_MIN;
                constexpr int LimitMax = EdgeInequalityIterator<TWarp, SubpixelRes>::LIMIT_MAX;

                int min_x = LimitMax;
                int max_x = LimitMin;
                for (int i = 0; i < TILE_SIZE / TWarp::SIZE; ++i) {
                    typename TWarp::Int
                        xmin,
                        xmax,
                        xmin0,
                        xmax0,
                        xmin1,
                        xmax1,
                        xmin2,
                        xmax2;

                    edge0.next_intervals(xmin0, xmax0);
                    edge1.next_intervals(xmin1, xmax1);
                    edge2.next_intervals(xmin2, xmax2);

                    //Intersection
                    
                    xmin0 = TWarp::max(xmin0, xmin1);
                    xmax0 = TWarp::min(xmax0, xmax1);

                    auto empty = TWarp::greater(xmin0, xmax0);
                    xmin0 = TWarp::select(empty, TWarp::scalar(LimitMax), xmin0);
                    xmax0 = TWarp::select(empty, TWarp::scalar(LimitMin), xmax0);

                    xmin = TWarp::max(xmin0, xmin2);
                    xmax = TWarp::min(xmax0, xmax2);

                    empty = TWarp::greater(xmin, xmax);
                    xmin = TWarp::select(empty, TWarp::scalar(LimitMax), xmin);
                    xmax = TWarp::select(empty, TWarp::scalar(LimitMin), xmax);

                    //To exclusive
                    xmax += 1;

                    TWarp::store(xmin, left + (b * TILE_SIZE + TWarp::SIZE * i));
                    TWarp::store(xmax, right + (b * TILE_SIZE + TWarp::SIZE * i));
                    
                    int l = TWarp::reduce_min(xmin);
                    int r = TWarp::reduce_max(xmax);
                    
                    min_x = l < min_x ? l : min_x;
                    max_x = r > max_x ? r : max_x;
                }
                t_left[b] = div_floor(min_x, TILE_SIZE);
                t_right[b] = div_ceil(max_x, TILE_SIZE);
                ++b;
                if (ty + b >= ty_end)
                    break;
            }

            auto column_step = TWarp::scalar(-1.f * TILE_SIZE / TWarp::SIZE);
            auto warp_size = TWarp::scalar(TWarp::SIZE);

            for (int k = 0; k < b; ++k, ++ty) {
                int tx_start = t_left[k];
                int tx_end = t_right[k];
                auto px_left = left + (k * TILE_SIZE);
                auto px_right = right + (k * TILE_SIZE);

                int sy = ty * TILE_SIZE;

                alignas(TWarp::ALIGNMENT) typename TWarp::Float attrs[NUM_ATTRS];
           
                for (int tx = tx_start; tx < tx_end; ++tx) {
                    int sx = tx * TILE_SIZE;

                    for (int i = 0; i < NUM_ATTRS; ++i)
                        attrs[i] = TWarp::fma(TWarp::scalar(sy * 1.f), a_dy[i], TWarp::fma(TWarp::scalar(sx * 1.f), a_dx[i], a_start[i]));

                    tile_lock(tx, ty);
                    for (int row = 0; row < TILE_SIZE; ++row) {
                        auto l = TWarp::scalar(px_left[row]);
                        auto r = TWarp::scalar(px_right[row]);
                        auto x = TWarp::add(TWarp::scalar(sx), idx);

                        for (int column = 0; column < TILE_SIZE; column += TWarp::SIZE) {
                            auto mask = TWarp::andnot(TWarp::greater_equals(x, l), TWarp::greater_equals(x, r));

                            fragment(tx, ty, column, row, mask, attrs);

                            for (int i = 0; i < NUM_ATTRS; ++i)
                                attrs[i] = TWarp::add(attrs[i], a_dx_wide[i]);

                            x = TWarp::add(x, warp_size);
                        }

                        for (int i = 0; i < NUM_ATTRS; ++i)
                            attrs[i] = TWarp::add(attrs[i], TWarp::fma(a_dx_wide[i], column_step, a_dy[i]));
                    }
                    tile_unlock(tx, ty);

                }
            }
        }
    }
}