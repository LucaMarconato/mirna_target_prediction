#ifndef TIMER_H
#define TIMER_H

#include <ctime>
#include <iostream>

// for benchmarking, this class is tremendously not thread-safe

extern clock_t my_time;

class Timer
{
public:
    static void start()
    {
        my_time = clock();
    }
    static void stop()
    {
        double t=static_cast<double>(clock()-my_time)/CLOCKS_PER_SEC;
        std::cerr << "elapsed: t = " << t << "\n";
    }
};

#endif //TIMER_H
