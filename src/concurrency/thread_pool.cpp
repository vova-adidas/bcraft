#include "thread_pool.hpp"

ThreadPool::ThreadPool(size_t numThreads) {
    start(numThreads);
}

ThreadPool::ThreadPool(): ThreadPool(std::thread::hardware_concurrency()) {
}

ThreadPool::~ThreadPool() {
    stop();
}

void ThreadPool::start(size_t numThreads) {
    for (size_t i = 0; i < numThreads; ++i) {
        m_workers.emplace_back([this] {
            Task task;
            while (m_tasks.consume_sync(task))
                task();
        });
    }
}

void ThreadPool::stop() {
    m_tasks.finish();

    for (auto& t : m_workers) {
        if (t.joinable())
            t.join();
    }
}
