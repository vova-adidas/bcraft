#pragma once
#include <vector>
#include <mutex>
#include <condition_variable>
#include <memory>
#include <cstddef>
#include <algorithm>

class MemoryPool {
public:
    MemoryPool(size_t total_size, size_t block_size)
        : m_total_size(total_size)
        , m_block_size(block_size) {

        m_block_count = total_size / block_size;
        m_buffer = std::unique_ptr<char[]>(new char[total_size]);
        m_free_blocks.reserve(m_block_count);
        for (size_t i = 0; i < m_block_count; ++i) {
            m_free_blocks.push_back(m_buffer.get() + i * block_size);
        }
    }

    std::shared_ptr<void> alloc() {
        std::unique_lock lock(m_mutex);
        m_cv.wait(lock, [this] { return !m_free_blocks.empty(); });

        void* ptr = m_free_blocks.back();
        m_free_blocks.pop_back();

        return {
            ptr,
            [this](void* p) { this->dealloc(p); }
        };
    }

    constexpr size_t block_size() const noexcept { return m_block_size; }
    size_t block_count() const noexcept { return m_block_count; }

private:

    void dealloc(void* p) {
        std::unique_lock lock(m_mutex);
        m_free_blocks.push_back(p);
        lock.unlock();
        m_cv.notify_one();
    }

    size_t m_total_size;
    size_t m_block_size;
    size_t m_block_count;

    std::unique_ptr<char[]> m_buffer;
    std::vector<void*> m_free_blocks;

    std::mutex m_mutex;
    std::condition_variable m_cv;
};