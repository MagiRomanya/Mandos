#ifndef MEMORY_POOL_H_
#define MEMORY_POOL_H_

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <new>

class MemoryPool {
    public:
        MemoryPool(size_t capacity) : capacity(capacity) {
            memory = malloc(capacity);
            if (!memory) {
                std::cerr << "ERROR::MemoryPool: Failed to allocate memory pool" << std::endl;
            }
            size = 0;
        }
        ~MemoryPool() {
            free(memory);
        }

        inline void reset() { size = 0; }

        inline void* allocate(size_t bytes) {
#ifndef NDEBUG
            if (size + bytes > capacity) {
                std::cerr << "ERROR::MemoryPool: out of memory!" << std::endl;
                throw std::bad_alloc();
            }
#endif
            void* ptr = (char*)(memory) + size;
            size += bytes;
            return ptr;
        }

    private:
        size_t capacity;
        size_t size = 0;
        void* memory;
};

template <typename T>
class MemoryPoolAllocator {
    public:
        MemoryPoolAllocator(MemoryPool& pool) : pool(pool) {}

        using value_type = T;

        inline T* allocate(std::size_t n) {
            return static_cast<T*>(pool.allocate(n * sizeof(T)));
        }

        inline void deallocate(T* p, size_t n) {
            // The pool manages it's own memory
        }

    private:
        MemoryPool& pool;
};

#endif // MEMORY_POOL_H_
