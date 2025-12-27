#pragma once

#include <mutex>
#include <condition_variable>

#include <queue>
#include <utility>

template<class T>
class ConcurrentQueue {
public:
	ConcurrentQueue() = default;

	~ConcurrentQueue() {
		finish();
	}

	void produce(T&& item) {

		std::lock_guard lock(m_mutex);

		m_queue.push(std::move(item));
		m_cond_var.notify_one();

	}

	auto size() {
		std::lock_guard lock(m_mutex);

		return m_queue.size();
	}

	bool consume(T& item) {

		std::lock_guard lock(m_mutex);

		if (m_queue.empty()) {
			return false;
		}

		item = std::move(m_queue.front());
		m_queue.pop();

		return true;

	}

	bool consume_sync(T& item) {

		std::unique_lock lock(m_mutex);

		m_sync_counter++;

		m_cond_var.wait(lock, [&] {
			return !m_queue.empty() || m_finish_processing;
		});

		if (m_queue.empty()) {
			if (--m_sync_counter == 0)
				m_sync_wait.notify_one();
			return false;
		}

		item = std::move(m_queue.front());
		m_queue.pop();

		if (--m_sync_counter == 0)
			m_sync_wait.notify_one();
		return true;

	}

	void finish() {

		std::unique_lock lock(m_mutex);

		m_finish_processing = true;
		m_cond_var.notify_all();

		m_sync_wait.wait(lock, [&]() {
			return m_sync_counter == 0;
		});

		m_finish_processing = false;

	}
private:

	std::queue<T> m_queue;

	std::mutex m_mutex;
	std::condition_variable m_cond_var;

	std::condition_variable m_sync_wait;
	bool m_finish_processing = false;
	int m_sync_counter = 0;
};