#include "chunk_generation.hpp"

#include <FastNoiseLite.h>

void chunk_generate_terrain(Chunk &chunk) {
    FastNoiseLite noise{};
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);

    int world_x = chunk.x() * Chunk_width;
    int world_z = chunk.z() * Chunk_width;

    auto blocks = chunk.blocks();

    for (int k = 0; k < Chunk_width; ++k) {

        if (chunk.is_destroyed())
            continue;

        for (int i = 0; i < Chunk_width; ++i) {
            auto x = world_x + i;
            auto z = world_z + k;
            int b = 0;
            for (int j = 64 + noise.GetNoise(x * 1.f, z * 1.f) * 5; j >= 0; j--) {
                auto y = j;

                if (b == 0)
                    blocks[i, j, k] = noise.GetNoise((float)x * 1, (float)y * 1, (float)z * 1) > 0 ? Blocks::Grass : Blocks::Void;
                else if (b < 4)
                    blocks[i, j, k] = noise.GetNoise((float)x * 1, (float)y * 1, (float)z * 1) > 0 ? Blocks::Dirt : Blocks::Void;
                else
                    blocks[i, j, k] = noise.GetNoise((float)x * 1, (float)y * 1, (float)z * 1) > 0 ? Blocks::Stone : Blocks::Void;

                b++;
            }
        }
    }
}
