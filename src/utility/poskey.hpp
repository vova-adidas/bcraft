#pragma once

struct PosKey {
	int x, y, z;

	bool operator==(const PosKey& other) const {
		return x == other.x && y == other.y && z == other.z;
	}
};

namespace std {
	template <>
	struct hash<PosKey> {
		std::size_t operator()(const PosKey& k) const noexcept {
			std::size_t h1 = std::hash<int>()(k.x);
			std::size_t h2 = std::hash<int>()(k.y);
			std::size_t h3 = std::hash<int>()(k.z);

			std::size_t seed = h1;
			seed ^= h2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			seed ^= h3 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			return seed;
		}
	};
}