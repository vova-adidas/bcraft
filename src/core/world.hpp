#pragma once 

#include "../concurrency/thread_pool.hpp"
#include "../utility/poskey.hpp"
#include "actor.hpp"

#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <utility>

struct Chunk;


class World {
public:

	constexpr static int Scan_size = 25;
	constexpr static int Max_actors = 10;

	World(int load_distance);

	bool add_actor(Actor* actor);

	void remove_actor(Actor* actor);

	void tick();

	int load_distance() const;

	int unload_distance() const;
private:

	ThreadPool thread_pool{};
	std::unordered_set<Actor*> m_actors{};
	std::unordered_map<Actor*, std::unordered_map<PosKey, bool>> m_occupied_chunks{};
	std::vector<std::pair<double, std::shared_ptr<Chunk>>> m_sorted_chunks{};
	std::unordered_map<PosKey, std::shared_ptr<Chunk>> m_chunks{};
	std::unordered_set<PosKey> m_occupied_positions{};
	int m_load_distance;
};