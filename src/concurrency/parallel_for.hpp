#pragma once

#include <algorithm>
#include <iterator>

#ifdef BCRAFT_STDEXECUTION
    #include <execution>
#elif defined(BCRAFT_OPENMP)
    #include <omp.h>
#endif

inline void parallel_foreach(auto iter_begin, auto iter_end, auto func) {
    #ifdef BCRAFT_STDEXECUTION
        std::for_each(std::execution::par, iter_begin, iter_end, func);
    #elif BCRAFT_OPENMP

        auto n = std::distance(iter_begin, iter_end);
        using diff_t = decltype(n);

        #pragma omp parallel for schedule(dynamic, 4)
        for (diff_t i = 0; i < n; ++i)
            func(*(iter_begin + i));

    #else 
        std::for_each(iter_begin, iter_end, func);
    #endif

}