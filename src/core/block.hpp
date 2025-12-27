#pragma once

class Block {
public:
	constexpr bool is_transparent() const { return *this == 0; }

	constexpr bool is_visible() const { return *this != 0; }

	constexpr bool is_opaque() const { return !is_transparent(); }

	constexpr bool operator==(Block& o) const { return m_id == o.m_id; }

	constexpr bool operator!=(Block& o) const { return m_id != o.m_id; }

	constexpr operator int() const { return m_id; }

	constexpr Block(int id) : m_id(id) {}

	constexpr Block() : m_id(0) {}

private:

	int m_id;
};

class Blocks {
public:
	static constexpr Block Void = { 0 };
	static constexpr Block Stone = { 1 };
	static constexpr Block Dirt = { 2 };
	static constexpr Block Grass = { 3 };
};