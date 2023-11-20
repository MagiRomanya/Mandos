#include <chrono>

class Clock {
public:
    Clock(double& ms_time) : ms_time(ms_time) {
        clock_begin = std::chrono::high_resolution_clock::now();
    }

    ~Clock() {
        auto clock_end = std::chrono::high_resolution_clock::now();
        ms_time = std::chrono::duration<double>(clock_end - clock_begin).count() * 1000;
    }

private:
    double& ms_time;
    std::chrono::high_resolution_clock::time_point clock_begin;
};
