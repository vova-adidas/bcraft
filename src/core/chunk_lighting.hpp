#pragma once

#include <array>

#include "chunk.hpp"

inline constexpr int Max_light = 15;

void chunk_fill_sky_light(Chunk& chunk);

void chunk_propagate_light_local(Chunk& chunk);

void chunk_finalize_light_propagation(const std::array<Chunk*, 9> &neighborhood);