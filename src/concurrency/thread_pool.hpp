#pragma once

#include <vector>
#include <thread>
#include <functional>

#include "concurrent_queue.hpp"

class ThreadPool {
public:
    using Task = std::function<void()>;

    explicit ThreadPool(size_t numThreads);

    ThreadPool();

    ~ThreadPool();

    void submit(Task&& task) {
        m_tasks.produce(std::move(task));
    }   

private:
    ConcurrentQueue<Task> m_tasks;
    std::vector<std::thread> m_workers;

    void start(size_t numThreads);

    void stop();
};