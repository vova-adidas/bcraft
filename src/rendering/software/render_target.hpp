#pragma once 

#include <xsimd/xsimd.hpp>
#include <vector>
#include <algorithm>

#include "../../concurrency/spinlock.hpp"
#include "warp.hpp"
#include "render_constants.hpp"

class RenderTarget {
public:
	RenderTarget(int width, int height) {
		resize(width, height);
	}

	void resize(int new_width, int new_height) {
		m_width = new_width;
		m_height = new_height;

		m_buf_width = static_cast<int>(std::ceil((new_width + 1.0) / TILE_SIZE)) * TILE_SIZE;
		m_buf_height = static_cast<int>(std::ceil((new_height + 1.0) / TILE_SIZE)) * TILE_SIZE;

		std::size_t size = m_buf_width * m_buf_height;

		if (std::size(m_visbuf_primitive_id) == size)
			return;

		m_tiled_width = (m_buf_width / TILE_SIZE);

		int tiled_size = m_tiled_width * (m_buf_height / TILE_SIZE);
		
		m_visbuf_primitive_id.resize(size);
		m_visbuf_weight_0.resize(size);
		m_visbuf_weight_1.resize(size);
		m_depth.resize(size);
		m_depth_mipmap.resize(tiled_size);
		m_spinlocks = std::make_unique<Spinlock[]>(tiled_size);
	}

	Warp::Int read_primitive(int x, int y) const {
		return Warp::Int::load_aligned(&m_visbuf_primitive_id[x + y * m_buf_width]);
	}

	void write_primitive(int x, int y, Warp::Int value) {
		return value.store_aligned(&m_visbuf_primitive_id[x + y * m_buf_width]);
	}

	Warp::Float read_weight_0(int x, int y) const {
		return Warp::Float::load_aligned(&m_visbuf_weight_0[x + y * m_buf_width]);
	}

	void write_weight_0(int x, int y, Warp::Float value) {
		return value.store_aligned(&m_visbuf_weight_0[x + y * m_buf_width]);
	}

	Warp::Float read_weight_1(int x, int y) const {
		return Warp::Float::load_aligned(&m_visbuf_weight_1[x + y * m_buf_width]);
	}

	void write_weight_1(int x, int y, Warp::Float value) {
		return value.store_aligned(&m_visbuf_weight_1[x + y * m_buf_width]);
	}

	Warp::Float read_depth(int x, int y) const {
		return Warp::Float::load_aligned(&m_depth[x + y * m_buf_width]);
	}

	void write_depth(int x, int y, Warp::Float value) {
		return value.store_aligned(&m_depth[x + y * m_buf_width]);
	}

	float read_depth_mipmap(int tx, int ty) const {
		return m_depth_mipmap[tx + ty * m_tiled_width];
	}

	void write_depth_mipmap(int tx, int ty, float value) {
		m_depth_mipmap[tx + ty * m_tiled_width] = value;
	}

	void lock_tile(int tx, int ty) {
		m_spinlocks[tx + ty * m_tiled_width].lock();
	}

	void unlock_tile(int tx, int ty) {
		m_spinlocks[tx + ty * m_tiled_width].unlock();
	}

	void clear() {
		std::fill(std::begin(m_visbuf_primitive_id), std::end(m_visbuf_primitive_id), -1);
		std::fill(std::begin(m_depth), std::end(m_depth), 0.0f);
		std::fill(std::begin(m_depth_mipmap), std::end(m_depth_mipmap), 0.0f);
	}

	int width() const { return m_width; }

	int height() const { return m_height; }

private:
	int m_width{}, m_height{}, m_buf_width{}, m_buf_height{}, m_tiled_width{};

	std::vector<int, xsimd::aligned_allocator<int, Warp::ALIGNMENT>> m_visbuf_primitive_id{};
	std::vector<float, xsimd::aligned_allocator<float, Warp::ALIGNMENT>> m_visbuf_weight_0{};
	std::vector<float, xsimd::aligned_allocator<float, Warp::ALIGNMENT>> m_visbuf_weight_1{};
	std::vector<float, xsimd::aligned_allocator<float, Warp::ALIGNMENT>> m_depth{};
	std::vector<float, xsimd::aligned_allocator<float, Warp::ALIGNMENT>> m_depth_mipmap{};
	std::unique_ptr<Spinlock[]> m_spinlocks = nullptr;

};