#pragma once 

#include <execution>
#include <memory>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <future>

#include <glm/glm.hpp>
#include <bgl/bgl.hpp>
#include <utility>

#include "../utility/poskey.hpp"
#include "../core/chunk.hpp"

#include "../concurrency/thread_pool.hpp"
#include "meshing.hpp"

#include "mesh_buffer_pool.hpp"

struct Vertex {
    glm::vec3 pos;
    int texture;
    int ao;
    int lighting;
};

constexpr int Subchunks_count = Chunk_height / 16;
constexpr int Subchunk_height = Chunk_height / Subchunks_count;

struct MeshData {
    Vertex vertices[Chunk_width * Chunk_width * Subchunk_height * 4 * 6];
    int indices[Chunk_width * Chunk_width * Subchunk_height * 6 * 9];
    std::size_t num_vertices;
    std::size_t num_indices;
};

class Renderer {
public:
    virtual ~Renderer() = default;

    virtual void* new_render_unit(
        std::shared_ptr<MeshData> data,
        int x,
        int y,
        int z,
        int local_min_x,
        int local_min_y,
        int local_min_z,
        int local_max_x,
        int local_max_y,
        int local_max_z) = 0;

    virtual void delete_render_unit(void* unit) = 0;
};

inline int block_to_texture_id(Meshing::Face face, Block block) {
    switch (block) {
        case Blocks::Stone:
            return 1;
        case Blocks::Grass:
            switch (face) {
                case Meshing::Face::Top: return 0;
                case Meshing::Face::Bot: return 2;
                default: return 3;
            }
            
        case Blocks::Dirt:
            return 2;
        default:
            return 0;
    }

}

struct ChunkRenderHandle;

class BuildMeshTask {
public:
    std::array<std::shared_ptr<Chunk>, 9> neighborhood;

    std::atomic<bool> is_canceled{};
    
    std::array<void*, Subchunks_count> result {};

    std::shared_ptr<MemoryPool> memory_pool {};
    
    std::shared_ptr<Renderer> renderer {};

    BuildMeshTask() = default;

    BuildMeshTask(std::shared_ptr<Renderer> renderer, std::shared_ptr<MemoryPool> memory_pool, std::array<std::shared_ptr<Chunk>, 9> neighborhood)
        : neighborhood(std::move(neighborhood)), is_canceled(false), memory_pool(std::move(memory_pool)), renderer(std::move(renderer)) {

    }

    void cancel() {
        is_canceled = true;
    }

    inline void run();
};

struct ChunkRenderHandle {

    std::shared_ptr<Renderer> renderer{};

    std::shared_ptr<BuildMeshTask> build_task{};

    std::shared_ptr<Chunk> chunk{};

    std::array<void*, Subchunks_count> render_units {};

    int neighbors_count = 0;

    std::atomic<bool> is_locked{};
    std::atomic<bool> is_dirty{};


    ~ChunkRenderHandle() {
        if (!renderer)
            return;

        for (int i = 0; i < Subchunks_count; ++i)
            if (render_units[i])
                renderer->delete_render_unit(render_units[i]);
    }
};

inline void BuildMeshTask::run() {

    int y_offset = 0;

    glm::vec3 v_min;
    glm::vec3 v_max;

    auto get_block = [&](int x, int y, int z) -> Block {

        y += y_offset;

        int dx = 1;
        int dz = 1;

        if (x < 0) {
            dx = 0;
            x = Chunk_width - 1;
        }
        if (x >= Chunk_width) {
            dx = 2;
            x = 0;
        }

        if (z < 0) {
            dz = 0;
            z = Chunk_width - 1;
        }
        if (z >= Chunk_width) {
            dz = 2;
            z = 0;
        }

        if (y < 0)
            return 0;
        if (y >= Chunk_height)
            return 0;

        return neighborhood[dx + dz * 3] != nullptr ? neighborhood[dx + dz * 3]->block_at(x, y, z) : Blocks::Void;
    };

    auto get_light = [&](int x, int y, int z) {
        y += y_offset;

        int dx = 1;
        int dz = 1;

        if (x < 0) {
            dx = 0;
            x = Chunk_width - 1;
        }
        if (x >= Chunk_width) {
            dx = 2;
            x = 0;
        }

        if (z < 0) {
            dz = 0;
            z = Chunk_width - 1;
        }
        if (z >= Chunk_width) {
            dz = 2;
            z = 0;
        }

        if (y < 0)
            return 0;
        if (y >= Chunk_height)
            return 15;

        return neighborhood[dx + dz * 3] != nullptr ? neighborhood[dx + dz * 3]->light_at(x, y, z) : 0;
    };

    Vertex* vertices;
    int* indices;
    int vdx;
    int idx;

    auto create_face = [&](auto f, auto face_dir, float* a, float* b, float* c, float* d) {
        const auto vertex_at = [&](int index) {
            if (index == 0)
                return glm::vec3(a[0], a[1], a[2]);
            if (index == 1)
                return glm::vec3(b[0], b[1], b[2]);
            if (index == 2)
                return glm::vec3(c[0], c[1], c[2]);
            return glm::vec3(d[0], d[1], d[2]);
        };

        int texture_id = block_to_texture_id(face_dir, f.block);
  
        int i = vdx;

        vertices[vdx++] = Vertex {
            vertex_at(0),
            texture_id,
            f.ao[0],
            f.light[0]
        };

        vertices[vdx++] = Vertex {
            vertex_at(1),
            texture_id,
            f.ao[1],
            f.light[1]
        };

        vertices[vdx++] = Vertex {
            vertex_at(2),
            texture_id,
            f.ao[2],
            f.light[2]
        };

        vertices[vdx++] = Vertex {
            vertex_at(3),
            texture_id,
            f.ao[3],
            f.light[3]
        };

        v_min = glm::min(v_min, vertex_at(0));
        v_min = glm::min(v_min, vertex_at(1));
        v_min = glm::min(v_min, vertex_at(2));
        v_min = glm::min(v_min, vertex_at(3));

        v_max = glm::max(v_max, vertex_at(0));
        v_max = glm::max(v_max, vertex_at(1));
        v_max = glm::max(v_max, vertex_at(2));
        v_max = glm::max(v_max, vertex_at(3));

        bool ao_flip = f.ao[0] + f.ao[3] > f.ao[1] + f.ao[2];
        if (ao_flip) {
            indices[idx++] = i;
            indices[idx++] = i + 1;
            indices[idx++] = i + 2;

            indices[idx++] = i + 2;
            indices[idx++] = i + 1;
            indices[idx++] = i + 3;
        }
        else {
            indices[idx++] = i;
            indices[idx++] = i + 1;
            indices[idx++] = i + 3;

            indices[idx++] = i + 3;
            indices[idx++] = i + 2;
            indices[idx++] = i + 0;
        }
    };

    auto& chunk = *neighborhood[4];

    for (int i = 0; i < Subchunks_count; ++i) {

        if (is_canceled)
            break;

        v_min = glm::vec3(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        v_max = -v_min;

        auto mesh_data = std::static_pointer_cast<MeshData>(memory_pool->alloc());

        vertices = mesh_data->vertices;
        indices = mesh_data->indices;
        idx = 0;
        vdx = 0;
        mesh_data->num_vertices = 0;
        mesh_data->num_indices = 0;
        Meshing::chunk_mesh<Chunk_width, Subchunk_height, Chunk_width>(get_block, get_light, create_face, is_canceled);

        y_offset += Subchunk_height;

        if (vdx == 0)
            continue;        

        if (is_canceled)
            break;

        mesh_data->num_vertices = vdx;
        mesh_data->num_indices = idx;

        result[i] = renderer->new_render_unit(
            mesh_data, 
            chunk.x() * Chunk_width, 
            chunk.y() + i * Subchunk_height,
            chunk.z() * Chunk_width, 
            v_min.x, v_min.y, v_min.z, 
            v_max.x, v_max.y, v_max.z);
    }
}


class View {
public:

    constexpr static int Block_size = sizeof(MeshData);

    explicit View(std::shared_ptr<Renderer> renderer) :
        m_framebuffer(nullptr), m_memory_pool(std::make_shared<MemoryPool>(sizeof(MeshData) * std::thread::hardware_concurrency() * 4, Block_size)), m_renderer(std::move(renderer)) {}

    void set_render_target(int* framebuffer, int width, int height) {
        m_framebuffer = framebuffer;
    }

    void push_chunk(const std::shared_ptr<Chunk>& chunk) {
        auto handle = std::make_shared<ChunkRenderHandle>();
        handle->chunk = chunk;
        handle->renderer = m_renderer;
        m_chunks.emplace(PosKey{ chunk->x(), chunk->y(), chunk->z() }, std::move(handle));
	}

    void update() {
        
        std::erase_if(m_chunks, [](auto& it) { return it.second->chunk->is_destroyed(); });
     
        for (auto& [pos, handle] : m_chunks) {
 
            auto& chunk = handle->chunk;

            std::array<std::shared_ptr<Chunk>, 9> neighborhood;
            int neighborhood_size = 0;

            for (int nx = -1; nx <= 1; ++nx)
                for (int nz = -1; nz <= 1; ++nz) {
                    auto neighbor = m_chunks.find({ chunk->x() + nx, 0, chunk->z() + nz });
                    if (neighbor != std::end(m_chunks)) {
                        neighborhood[nx + 1 + (nz + 1) * 3] = neighbor->second->chunk;
                        neighborhood_size++;
                    }
                }

            if (handle->neighbors_count < neighborhood_size)
                handle->is_dirty = true;

            handle->neighbors_count = neighborhood_size;


            if (handle->is_dirty) {

                if (handle->build_task != nullptr) {
                    handle->build_task->cancel();
                    handle->build_task = nullptr;
                }

                if (handle->is_locked)
                    continue;

                handle->build_task = std::make_shared<BuildMeshTask>(m_renderer, m_memory_pool, neighborhood);
                handle->is_locked = true;
                handle->is_dirty = false;

                thread_pool.submit([handle = handle, task = handle->build_task] {
                    task->run();

                    for (int i = 0; i < Subchunks_count; ++i) {
                        //move

                        if (handle->render_units[i] != nullptr)
                            handle->renderer->delete_render_unit(handle->render_units[i]);

                        handle->render_units[i] = task->result[i];
                    }

                    handle->is_locked = false;
                });

            }
        }
    }

private:
    
    std::shared_ptr<Renderer> m_renderer{};

    ThreadPool thread_pool{};

    std::shared_ptr<MemoryPool> m_memory_pool;

    std::unordered_map<PosKey, std::shared_ptr<ChunkRenderHandle>> m_chunks{};
    
    int* m_framebuffer = nullptr;
};