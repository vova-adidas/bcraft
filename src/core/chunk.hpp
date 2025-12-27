#pragma once

#include "block.hpp"

#include <atomic>

constexpr int Chunk_width = 16;
constexpr int Chunk_height = 256;

enum class ChunkState {
	Empty,
	LocalReady,
	Ready
};

struct Chunk {

	constexpr static int Buffer_Size = Chunk_width * Chunk_width * Chunk_height;

	template <typename T>
	class View {
	public:
		explicit View(T data) : m_data(data) {}

		const auto& operator[](int x, int y, int z) const { return m_data[idx(x, y, z)]; }
		auto& operator[](int x, int y, int z) { return m_data[idx(x, y, z)]; }

		auto* data() { return m_data; }
		const auto* begin() const { return m_data; }
		const auto* end() const { return m_data + Buffer_Size; }
		auto* begin() { return m_data; }
		auto* end() { return m_data + Buffer_Size; }

	private:
		T m_data;
	};

	Chunk(int x, int y, int z) : m_x(x), m_y(y), m_z(z) {}

	Chunk(const Chunk&) = delete;

	Chunk& operator=(const Chunk&) = delete;

	Block& block_at(int x, int y, int z) { return m_blocks[idx(x, y, z)]; }
	Block block_at(int x, int y, int z) const { return m_blocks[idx(x, y, z)]; }

	int& light_at(int x, int y, int z) { return m_lighting[idx(x, y, z)]; }
	int light_at(int x, int y, int z) const { return m_lighting[idx(x, y, z)]; }

	void mark_destroyed() { m_is_destroyed = true; }
	bool is_destroyed() const { return m_is_destroyed; }

	void set_state(ChunkState state) { m_state = state; }
	ChunkState state() const { return m_state; }

	int& x() { return m_x; }
	int x() const { return m_x; }

	int& y() { return m_y; }
	int y() const { return m_y; }

	int& z() { return m_z; }
	int z() const { return m_z; }

	View<Block*> blocks() { return View { m_blocks }; }

	View<const Block*> blocks() const { return View { m_blocks }; }

	View<int*> lighting() { return View { m_lighting }; }

	View<const int*> lighting() const { return View {m_lighting}; }

	void lock() { m_is_locked = true; }

	void unlock() { m_is_locked = false; }

	bool is_locked() const { return m_is_locked; }

private:
	static int idx(int x, int y, int z) { return x + Chunk_width * (y + z * Chunk_height); }

	Block m_blocks[Buffer_Size] {};
	int m_lighting[Buffer_Size] {};

	int m_x, m_y, m_z;
	
	std::atomic<ChunkState> m_state {};
	std::atomic<bool> m_is_destroyed {};
	std::atomic<bool> m_is_locked {};
};