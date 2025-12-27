#pragma once

class Meshing {
public:

    template<int Width, int Height>
    static void greedy(auto is_visible, auto get_face_data, auto set_face) {
        using FaceData = decltype(get_face_data(0, 0));

        int w[Width];
        int c[Width];
        FaceData fds[Width];

        const auto row_exist_at = [&](int x) {
            return w[x] != 0;
        };

        const auto place_row_at = [&](int x, int y) {
            set_face(fds[x], x, y - c[x], w[x], c[x]);
            w[x] = 0;
        };

        const auto place_row_at2 = [&](int x, int y, int& row_length, const FaceData& previous_fd) {
            if (row_length == 0)
                return;

            auto row_start = x - row_length;
            if (w[row_start] == row_length && fds[row_start] == previous_fd)
                c[row_start]++;
            else
            {
                if (row_exist_at(row_start))
                    place_row_at(row_start, y);

                w[row_start] = row_length;
                c[row_start] = 1;
                fds[row_start] = previous_fd;
            }
            row_length = 0;
        };

        for (int i = 0; i < Width; ++i)
            w[i] = 0;

        for (int y = 0; y < Height; ++y) {
            int row_length = 0;

            alignas(FaceData) char b[sizeof(FaceData)];
            auto* previous_fd = (FaceData*)&b;

            for (int x = 0; x < Width; ++x) {

                if (is_visible(x, y)) {
                    auto fd = get_face_data(x, y);

                    if (!(*previous_fd == fd))
                        place_row_at2(x, y, row_length, *previous_fd);

                    if (row_length == 0)
                        row_length = 1;
                    else {
                        ++row_length;
                        if (row_exist_at(x))
                            place_row_at(x, y);
                    }
                    *previous_fd = fd;
                }
                else {
                    place_row_at2(x, y, row_length, *previous_fd);

                    if (row_exist_at(x))
                        place_row_at(x, y);
                }
            }

            place_row_at2(Width, y, row_length, *previous_fd);
        }

        for (int x = 0; x < Width; ++x) {
            int y = Height;
            if (w[x] != 0)
                set_face(fds[x], x, y - c[x], w[x], c[x]);
        }
    }

    enum class Face {
        Front,
        Back,
        Left,
        Right,
        Top,
        Bot
    };

    template<int Width, int Height, int Length>
    static void chunk_mesh(auto block_selector, auto light_selector, auto set_face, auto& is_canceled) {

        using Block = decltype(block_selector(0, 0, 0));

        struct FaceData {
            Block block;
            int ao[4];
            int light[4];

            bool operator==(const FaceData& other) const {

                bool cmp = block == other.block;
                if (!cmp)
                    return false;
                for (int i = 0; i < 4; ++i) {
                    cmp &= light[i] == other.light[i];
                    if (!cmp)
                        return false;
                }
                for (int i = 0; i < 4; ++i) {
                    cmp &= ao[i] == other.ao[i];
                    if (!cmp)
                        return false;
                }
                return cmp;
            }
        };

        auto compute_ao = [&](auto block_selector, int x, int y) {

            auto b0 = block_selector(x, y);
            auto b1 = block_selector(x, y - 1);
            auto b2 = block_selector(x - 1, y);
            auto b3 = block_selector(x - 1, y - 1);

            return (b0 != 0) + (b1 != 0) + (b2 != 0) + (b3 != 0);
        };

        auto compute_light = [&](auto light_selector, auto block_selector, int x, int y) -> int {

            std::pair<int, int> uv[] = {
                { x, y },
                { x, y - 1 },
                { x - 1, y },
                { x - 1, y - 1 },
            };

            auto total = 0;
            auto sum = 0;
            for (int i = 0; i < 4; ++i)
                if (block_selector(uv[i].first, uv[i].second).is_transparent()) {
                    sum += light_selector(uv[i].first, uv[i].second);
                    ++total;
                }

            if (total > 0)
                return static_cast<float>(sum) / total;

            return 0;

            };

        auto face_data = [&](bool ccw, auto block_selector, auto next_block_selector, auto light_selector, int x, int y) {
            FaceData fd;
            fd.block = block_selector(x, y);
            fd.light[0] = compute_light(light_selector, next_block_selector, x, y);
            fd.ao[0] = compute_ao(next_block_selector, x, y);
            if (ccw) {
                fd.light[1] = compute_light(light_selector, next_block_selector, x + 1, y);
                fd.ao[1] = compute_ao(next_block_selector, x + 1, y);
                fd.light[2] = compute_light(light_selector, next_block_selector, x, y + 1);
                fd.ao[2] = compute_ao(next_block_selector, x, y + 1);
            }
            else {
                fd.light[1] = compute_light(light_selector, next_block_selector, x, y + 1);
                fd.ao[1] = compute_ao(next_block_selector, x, y + 1);
                fd.light[2] = compute_light(light_selector, next_block_selector, x + 1, y);
                fd.ao[2] = compute_ao(next_block_selector, x + 1, y);
            }
            fd.light[3] = compute_light(light_selector, next_block_selector, x + 1, y + 1);
            fd.ao[3] = compute_ao(next_block_selector, x + 1, y + 1);
            //191.852, 11.6534, -455.004
            //for (int i = 0; i < 4; ++i)
            //    fd.light[i] = light_selector(x, y);
            return fd;
        };

        auto is_visible = [&](int x, int y, int z, int nx, int ny, int nz) {
            return block_selector(x, y, z).is_visible() && block_selector(nx, ny, nz).is_transparent();
        };

        float epsilon = 0.0001f;

        for (int z = 0; z < Length; ++z) {
            if (is_canceled)
                return;

            greedy<Width, Height>(
                [&](int x, int y) { return is_visible(x, y, z, x, y, z + 1); },
                [&](int x, int y) { return face_data(false, [&](int x, int y) { return block_selector(x, y, z); }, [&](int x, int y) { return block_selector(x, y, z + 1); }, [&](int x, int y) { return light_selector(x, y, z + 1); }, x, y); },
                [&](auto f, float x, float y, float w, float h) {
                    float a[]{ x - epsilon, y - epsilon, z + 1.f };
                    float b[]{ x - epsilon, y + h + epsilon, z + 1.f };
                    float c[]{ x + w + epsilon, y - epsilon, z + 1.f };
                    float d[]{ x + w + epsilon, y + h + epsilon, z + 1.f };

                    set_face(f, Face::Front, a, b, c, d);
            });

            greedy<Width, Height>(
                [&](int x, int y) { return is_visible(x, y, z, x, y, z - 1); },
                [&](int x, int y) { return face_data(true, [&](int x, int y) { return block_selector(x, y, z); }, [&](int x, int y) { return block_selector(x, y, z - 1); }, [&](int x, int y) { return light_selector(x, y, z - 1); }, x, y); },
                [&](auto f, float x, float y, float w, float h) {
                    float a[]{ x - epsilon, y - epsilon, z * 1.f };
                    float b[]{ x + w + epsilon, y - epsilon, z * 1.f };
                    float c[]{ x - epsilon, y + h + epsilon, z * 1.f };
                    float d[]{ x + w + epsilon, y + h + epsilon, z * 1.f };

                    set_face(f, Face::Back, a, b, c, d);
             });
        }

        for (int y = 0; y < Height; ++y) {
            if (is_canceled)
                return;

            greedy<Width, Length>(
                [&](int x, int z) { return is_visible(x, y, z, x, y + 1, z); },
                [&](int x, int z) { return face_data(true, [&](int x, int z) { return block_selector(x, y, z); }, [&](int x, int z) { return block_selector(x, y + 1, z); }, [&](int x, int z) { return light_selector(x, y + 1, z); }, x, z); },
                [&](auto f, float x, float z, float w, float h) {
                    float a[]{ x - epsilon, y + 1.f, z - epsilon };
                    float b[]{ x + w + epsilon, y + 1.f, z - epsilon };
                    float c[]{ x - epsilon, y + 1.f, z + h + epsilon };
                    float d[]{ x + w + epsilon, y + 1.f, z + h + epsilon };

                    set_face(f, Face::Top, a, b, c, d);
            });

            greedy<Width, Length>(
                [&](int x, int z) { return is_visible(x, y, z, x, y - 1, z);; },
                [&](int x, int z) { return face_data(false, [&](int x, int z) { return block_selector(x, y, z); }, [&](int x, int z) { return block_selector(x, y - 1, z); }, [&](int x, int z) { return light_selector(x, y - 1, z); }, x, z); },
                [&](auto f, float x, float z, float w, float h) {
                    float a[]{ x - epsilon, y * 1.f, z - epsilon };
                    float b[]{ x - epsilon, y * 1.f, z + h + epsilon };
                    float c[]{ x + w + epsilon, y * 1.f, z - epsilon };
                    float d[]{ x + w + epsilon, y * 1.f, z + h + epsilon };

                    set_face(f, Face::Bot, a, b, c, d);
            });
        }

        for (int x = 0; x < Width; ++x) {
            if (is_canceled)
                return;

            greedy<Height, Length>(
                [&](int y, int z) { return is_visible(x, y, z, x + 1, y, z); },
                [&](int y, int z) { return face_data(false, [&](int y, int z) { return block_selector(x, y, z); }, [&](int y, int z) { return block_selector(x + 1, y, z); }, [&](int y, int z) { return light_selector(x + 1, y, z); }, y, z); },
                [&](auto f, float y, float z, float w, float h) {
                    float a[]{ x + 1.f, y - epsilon, z};
                    float b[]{ x + 1.f, y - epsilon, z + h + epsilon };
                    float c[]{ x + 1.f, y + w + epsilon, z - epsilon };
                    float d[]{ x + 1.f, y + w + epsilon, z + h + epsilon };

                    set_face(f, Face::Right, a, b, c, d);
            });

            greedy<Height, Length>(
                [&](int y, int z) { return is_visible(x, y, z, x - 1, y, z); },
                [&](int y, int z) { return face_data(true, [&](int y, int z) { return block_selector(x, y, z); }, [&](int y, int z) { return block_selector(x - 1, y, z); }, [&](int y, int z) { return light_selector(x - 1, y, z); }, y, z); },
                [&](auto f, float y, float z, float w, float h) {
                    float a[]{ x * 1.f, y - epsilon, z - epsilon };
                    float b[]{ x * 1.f, y + w + epsilon, z - epsilon };
                    float c[]{ x * 1.f, y - epsilon, z + h + epsilon };
                    float d[]{ x * 1.f, y + w + epsilon, z + h + epsilon };

                    set_face(f, Face::Left, a, b, c, d);
            });
        }
    }
};