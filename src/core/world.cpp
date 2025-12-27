#include "world.hpp"
#include "chunk.hpp"
#include "chunk_lighting.hpp"
#include "chunk_generation.hpp"

#include <algorithm>

#include "bgl/raster.hpp"

World::World(int load_distance) : m_load_distance(load_distance) {
	m_occupied_positions.reserve((Scan_size + 1) * (Scan_size + 1) * Max_actors);
	m_sorted_chunks.reserve((Scan_size + 1) * (Scan_size + 1) * Max_actors);
}

bool World::add_actor(Actor *actor) {
    if (std::size(m_actors) < Max_actors) {
        m_actors.emplace(actor);
        m_occupied_chunks.emplace(actor, std::unordered_map<PosKey, bool>{});
        return true;
    }
    return false;
}

void World::remove_actor(Actor *actor) {
    m_actors.erase(actor);
    m_occupied_chunks.erase(actor);
}

void World::tick() {
	m_occupied_positions.clear();

	for (auto* actor : m_actors) {
		double actor_x = actor->pos.x;
		double actor_z = actor->pos.z;
		auto center_x = actor_x / Chunk_width;
		auto center_z = actor_z / Chunk_width;
		auto radius = static_cast<double>(load_distance());
		bgl::EllipseScanlineIteratorBase it(center_x, center_z, radius, radius);

		while (it.next()) {
			int z = it.line();
			// Мы берем floor/ceil от границ, чтобы захватить все задетые чанки
			int x_start = static_cast<int>(it.begin());
			int x_end   = static_cast<int>(it.end());

			for (int x = x_start; x <= x_end; ++x) {
				// Теперь нам не нужно считать мировые координаты и проверять дистанцию
				// Итератор уже гарантирует, что эти x и z внутри круга
				PosKey pos { x, 0, z };

				if (m_occupied_chunks[actor].find(pos) == m_occupied_chunks[actor].end()) {
					m_occupied_chunks[actor].emplace(pos, false);
				}
			}
		}

		std::erase_if(m_occupied_chunks[actor], [&](auto& item) {
			auto& [pos, _] = item;

			double chunk_world_center_x = pos.x * Chunk_width + Chunk_width / 2.;
			double chunk_world_center_z = pos.z * Chunk_width + Chunk_width / 2.;

			auto dir_x = actor_x - chunk_world_center_x;
			auto dir_z = actor_z - chunk_world_center_z;
			auto sq_distance = dir_x * dir_x + dir_z * dir_z;
			auto radius = unload_distance() * Chunk_width;
			return sq_distance > radius * radius;

		});

		for (auto& [pos, _] : m_occupied_chunks[actor])
			m_occupied_positions.emplace(pos);
	}

	for (auto pos : m_occupied_positions)
		if (m_chunks.find(pos) == std::end(m_chunks)) {
			auto ch = std::make_shared<Chunk>(pos.x, pos.y, pos.z);
			m_chunks.emplace(pos, std::move(ch));
		}

	std::erase_if(m_chunks, [&](auto& item) {
		auto& [pos, ch] = item;
		bool destroyed = m_occupied_positions.find(pos) == std::end(m_occupied_positions);
		if (destroyed)
			ch->mark_destroyed();
		return destroyed;
	});

	m_sorted_chunks.clear();
	for (auto& [pos, chunk] : m_chunks) {
		double chunk_world_center_x = pos.x * Chunk_width + Chunk_width / 2.;
		double chunk_world_center_z = pos.z * Chunk_width + Chunk_width / 2.;
		double min_distance = std::numeric_limits<double>::max();
		for (auto* actor : m_actors) {
			double actor_x = actor->pos.x;
			double actor_z = actor->pos.z;
			auto dir_x = actor_x - chunk_world_center_x;
			auto dir_z = actor_z - chunk_world_center_z;
			auto sq_distance = dir_x * dir_x + dir_z * dir_z;
			if (sq_distance < min_distance)
				min_distance = sq_distance;
		}

		m_sorted_chunks.emplace_back(min_distance, chunk);
	}

	std::sort(std::begin(m_sorted_chunks), std::end(m_sorted_chunks), [](auto& lhs, auto& rhs) { return lhs.first < rhs.first;});

	for (auto& [_, chunk] : m_sorted_chunks) {

		if (chunk->is_locked())
			continue;

		switch (chunk->state()) {
			case ChunkState::Empty: {
				chunk->lock();

				thread_pool.submit([chunk = chunk] {
					chunk_generate_terrain(*chunk);
					chunk_fill_sky_light(*chunk);
					chunk_propagate_light_local(*chunk);
					chunk->set_state(ChunkState::LocalReady);
					chunk->unlock();
				});


				break;
			}
			case ChunkState::LocalReady: {

				std::array<std::shared_ptr<Chunk>, 9> neighborhood;
				int neighborhood_size = 0;

				for (int nx = -1; nx <= 1; ++nx)
					for (int nz = -1; nz <= 1; ++nz) {
						auto neighbor = m_chunks.find({ chunk->x() + nx, 0, chunk->z() + nz });
						if (neighbor != std::end(m_chunks)) {
							neighborhood[nx + 1 + (nz + 1) * 3] = neighbor->second;
							if(neighbor->second->state() >= ChunkState::LocalReady)
								neighborhood_size++;
						}
					}

				if (neighborhood_size == 9) {
					chunk->lock();

					thread_pool.submit([chunk = chunk, neighborhood = neighborhood, neighborhood_size = neighborhood_size] {

						std::array<Chunk*, 9> neighborhood_raw {};
						for (int i = 0; i < neighborhood_size; ++i)
							neighborhood_raw[i] = neighborhood[i].get();

						chunk_finalize_light_propagation(neighborhood_raw);
						chunk->set_state(ChunkState::Ready);

						chunk->unlock();
					});
				}

				break;
			}
			default:
				break;
		}

	}

	for (auto&[pos, chunk] : m_chunks)
		if (chunk->state() == ChunkState::Ready)
			for (auto& [actor, occupied_chunks] : m_occupied_chunks) {
				auto it = occupied_chunks.find(pos);
				if (it != std::end(occupied_chunks) && !it->second) {
					//Запомнить, что уже отправлял
					it->second = true;
					actor->receive_chunk(chunk);
				}
			}

}

int World::load_distance() const {
	return m_load_distance;
}

int World::unload_distance() const {
	return load_distance() + 2;
}
