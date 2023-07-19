#include "WorkPool.h"


WorkPool::WorkPool(): last_idx_done(-1)
{
    pthread_mutex_init(&lock, NULL);
}

void WorkPool::new_partition()
{
    if (!works.size() || works.back().size())
        works.push_back(std::vector<uint64_t>());
}

void WorkPool::pop_back_empty_partition()
{
    if (!works.back().size())
        works.pop_back();
}

void WorkPool::push_work(uint64_t work_cell)
{
    works.back().push_back(work_cell);
}

const std::vector<uint64_t>* WorkPool::get_work()  //thread-safe
{
    uint64_t index_partition;
    pthread_mutex_lock(&lock);
    if (last_idx_done == works.size() - 1)
        index_partition = -1;
    else
        index_partition = ++last_idx_done;
    pthread_mutex_unlock(&lock);
    
    return (index_partition != -1) ? &(works.at(index_partition)) : nullptr;
}

void WorkPool::reset_work_completion()
{
    last_idx_done = -1;
}

void WorkPool::finalize()
{
    pthread_mutex_destroy(&lock);
}