#include "chunk_lighting.hpp"

#include <vector>

struct LightStep {
	int x, y, z, light;
};

class LightQueue {
public:
	constexpr static std::size_t Size = Chunk_width * Chunk_height * Chunk_width * 10;

	LightQueue() : m_queue(make_buffer()), m_head(0), m_tail(0) {
	
	}

	bool empty() const { return m_head == m_tail; }

	void push(int x, int y, int z, int light) {
		std::size_t next = (m_tail + 1) % Size;
		
		//if (next == head)
		//	return;
		
		m_queue[m_tail] = {x, y, z, light};
		m_tail = next;
	};

	LightStep pop() {
		auto s = m_queue[m_head];
		m_head = (m_head + 1) % Size;
		return s;
	};


private:
	static std::vector<LightStep>& make_buffer() {
		thread_local std::vector<LightStep> queue(Size);
		return queue;
	}

	std::vector<LightStep>& m_queue;
	std::size_t m_head, m_tail;
};


static void chunk_set_light(Chunk& chunk, int p_x, int p_y, int p_z, int new_light) {

	auto blocks = chunk.blocks();
	auto lightmap = chunk.lighting();

	if (blocks[p_x, p_y, p_z].is_opaque())
		return;

	if (new_light < lightmap[p_x, p_y, p_z])
		return;

	LightQueue queue{};

	queue.push(p_x, p_y, p_z, new_light);

	while (!queue.empty()) {
		auto [x, y, z, current_light] = queue.pop();
		lightmap[x, y, z] = current_light;
	
		current_light--;

		if (x + 1 < Chunk_width && current_light > lightmap[x + 1, y, z] && blocks[x + 1, y, z].is_transparent())
			queue.push(x + 1, y, z, current_light);
		if (x - 1 >= 0 && current_light > lightmap[x - 1, y, z] && blocks[x - 1, y, z].is_transparent())
			queue.push(x - 1, y, z, current_light);

		if (z + 1 < Chunk_width && current_light > lightmap[x, y, z + 1] && blocks[x, y, z + 1].is_transparent())
			queue.push(x, y, z + 1, current_light);
		if (z - 1 >= 0 && current_light > lightmap[x, y, z - 1] && blocks[x, y, z - 1].is_transparent())
			queue.push(x, y, z - 1, current_light);

		if (y + 1 < Chunk_height && current_light > lightmap[x, y + 1, z] && blocks[x, y + 1, z].is_transparent())
			queue.push(x, y + 1, z, current_light);
		if (y - 1 >= 0 && current_light > lightmap[x, y - 1, z] && blocks[x, y + 1, z].is_transparent())
			queue.push(x, y - 1, z, current_light);
	}

}

void chunk_fill_sky_light(Chunk& chunk) {
	auto blocks = chunk.blocks();
	auto lightmap = chunk.lighting();

	for (int x = 0; x < Chunk_width; ++x) {

		if (chunk.is_destroyed())
			continue;

		for (int z = 0; z < Chunk_width; ++z) {

			int light = Max_light;
			for (int y = Chunk_height - 1; y >= 0; --y) {
				lightmap[x, y, z] = light;

				if (blocks[x, y, z].is_opaque())
					light = 0;
			}
		}
	}
}

void chunk_propagate_light_local(Chunk& chunk) {
	auto lightmap = chunk.lighting();

	for (int x = 0; x < Chunk_width; ++x) {
		if (chunk.is_destroyed())
			continue;

		for (int z = 0; z < Chunk_width; ++z)
			for (int y = 0; y < Chunk_height; y++)
				if (lightmap[x, y, z] != 0)
					chunk_set_light(chunk, x, y, z, lightmap[x, y, z]);
	}
}

void chunk_finalize_light_propagation(const std::array<Chunk*, 9> &neighborhood) {

	const auto neighbor_at = [&](int x, int z) -> Chunk& { return *neighborhood[x + 1 + (z + 1) * 3]; };

	auto& chunk = neighbor_at(0, 0);


	auto neighbor_bot_lightmap = neighbor_at(0, -1).lighting();
	auto neighbor_top_lightmap = neighbor_at(0, 1).lighting();
	auto neighbor_bot_blocks = neighbor_at(0, -1).blocks();
	auto neighbor_top_blocks = neighbor_at(0, 1).blocks();

	for (int y = 0; y < Chunk_height; ++y)
		for (int x = 0; x < Chunk_width; ++x) {
			chunk_set_light(chunk, x, y, 0, neighbor_bot_lightmap[x, y, Chunk_width - 1] - 1);
			chunk_set_light(chunk, x, y, Chunk_width - 1, neighbor_top_lightmap[x, y, 0] - 1);
		}

	auto neighbor_left_lightmap = neighbor_at(-1, 0).lighting();
	auto neighbor_right_lightmap = neighbor_at(1, 0).lighting();
	auto neighbor_left_blocks = neighbor_at(-1, 0).blocks();
	auto neighbor_right_blocks = neighbor_at(1, 0).blocks();

	for (int y = 0; y < Chunk_height; ++y)
		for (int z = 0; z < Chunk_width; ++z) {
			chunk_set_light(chunk, 0, y, z, neighbor_left_lightmap[Chunk_width - 1, y, z] - 1);
			chunk_set_light(chunk, Chunk_width - 1, y, z, neighbor_right_lightmap[0, y, z] - 1);
		}

	auto neighbor_left_bot_lightmap = neighbor_at(-1, -1).lighting();

	for (int y = 0; y < Chunk_height; ++y) {
		auto l = neighbor_left_bot_lightmap[Chunk_width - 1, y, Chunk_width - 1];

		auto left_block = neighbor_left_blocks[Chunk_width - 1, y, 0];
		auto bottom_block = neighbor_bot_blocks[0, y, Chunk_width - 1];

		if (left_block.is_transparent() || bottom_block.is_transparent())
			chunk_set_light(chunk, 0, y, 0, l - 2);
	}

	auto neighbor_right_bot_lightmap = neighbor_at(1, -1).lighting();

	for (int y = 0; y < Chunk_height; ++y) {
		auto l = neighbor_right_bot_lightmap[0, y, Chunk_width - 1];

		auto right_block = neighbor_right_blocks[0, y, 0];
		auto bottom_block = neighbor_bot_blocks[Chunk_width - 1, y, Chunk_width - 1];

		if (right_block.is_transparent() || bottom_block.is_transparent())
			chunk_set_light(chunk, Chunk_width - 1, y, 0, l - 2);
	}

	auto neighbor_left_top_lightmap = neighbor_at(-1, 1).lighting();

	for (int y = 0; y < Chunk_height; ++y) {
		auto l = neighbor_left_top_lightmap[Chunk_width - 1, y, 0];

		auto left_block = neighbor_left_blocks[Chunk_width - 1, y, Chunk_width - 1];
		auto top_block = neighbor_top_blocks[0, y, 0];
		if (left_block.is_transparent() || top_block.is_transparent())
			chunk_set_light(chunk, 0, y, Chunk_width - 1, l - 2);
	}

	auto neighbor_right_top_lightmap = neighbor_at(1, 1).lighting();

	for (int y = 0; y < Chunk_height; ++y) {
		auto l = neighbor_right_top_lightmap[0, y, 0];

		auto right_block = neighbor_right_blocks[0, y, Chunk_width - 1];
		auto top_block = neighbor_top_blocks[Chunk_width - 1, y, 0];
		if (right_block.is_transparent() || top_block.is_transparent())
			chunk_set_light(chunk, Chunk_width - 1, y, Chunk_width - 1, l - 2);
	}
}