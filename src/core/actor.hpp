#pragma once

#include <memory>

struct Chunk;

struct Actor {
    struct {
        double x, y, z;
    } pos {};

    virtual void receive_chunk(std::shared_ptr<Chunk> event) = 0;

    virtual ~Actor() = default;
};
