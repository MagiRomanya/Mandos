#include "physics_state.hpp"
#include "simulation.hpp"
#include <iostream>
#include <thread>
#include <condition_variable>
#include <queue>

namespace
{

static std::mutex queue_mutex;
static std::mutex state_mutex;
static mandos::PhysicsState state_async;
static mandos::EnergyAndDerivatives f_async(0);
static std::condition_variable cv;
static std::queue<int> simulation_requests;
std::atomic<bool> stop_requested = false;
}  // namespace

namespace mandos
{

void __simulation_loop_async(const Simulation simulation)
{
    while (true) {
        // Check for simulation requests
        // ---------------------------------------------------------------
        std::unique_lock<std::mutex> queue_lock(queue_mutex);
        cv.wait(queue_lock, [] { return (!simulation_requests.empty()) or stop_requested; });
        // Stop thread condition
        if (stop_requested) {
            break;
        }

        int request = simulation_requests.front();
        simulation_requests.pop();
        queue_lock.unlock();

        // Preform the simulation step
        // ---------------------------------------------------------------
        std::lock_guard<std::mutex> state_lock(state_mutex);
        for (int i = 0; i < request; i++)
            simulation_step(simulation, state_async, f_async);
    }
}

/**
 * Struct that starts a thread on creation and stops and joins it on destruction.
 */
struct SimulationLoopThread {
    SimulationLoopThread(const Simulation simulation)
    {
        simulation_loop_thread = std::jthread(__simulation_loop_async, simulation);
    }
    ~SimulationLoopThread()
    {
        stop_requested = true;
        cv.notify_one();
        simulation_loop_thread.join();
    }

    std::jthread simulation_loop_thread;
};

static std::unique_ptr<SimulationLoopThread> simulation_loop_thread;

void simulation_async_loop(Simulation simulation)
{
    state_async = simulation.initial_state;
    simulation_loop_thread = std::make_unique<SimulationLoopThread>(simulation);
}

void simulation_async_loop_request_iteration()
{
    std::lock_guard<std::mutex> lock(queue_mutex);
    if (simulation_requests.size() == 0)
        simulation_requests.push(1);
    cv.notify_one();
}

void set_current_physics_state(PhysicsState state)
{
    std::lock_guard<std::mutex> lock(state_mutex);
    state_async = state;
}

PhysicsState get_current_physics_state()
{
    return state_async;
}
EnergyAndDerivatives get_current_energy_and_derivatives()
{
    return f_async;
}

}  // namespace mandos
