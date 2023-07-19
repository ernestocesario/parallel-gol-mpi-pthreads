#ifndef WORK_POOL_H
#define WORK_POOL_H

#include <stdint.h>
#include <pthread.h>
#include <vector>

class WorkPool
{
    private:
        std::vector<std::vector<uint64_t>> works;
        uint64_t last_idx_done;
        pthread_mutex_t lock;

    public:
        WorkPool();
        void new_partition();
        void pop_back_empty_partition();
        void push_work(uint64_t work_cell);
        const std::vector<uint64_t>* get_work();
        void reset_work_completion();
        void finalize();
};


#endif