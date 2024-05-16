#ifndef EA5CC9BC_9B26_421C_B4CB_4E265608D7F2
#define EA5CC9BC_9B26_421C_B4CB_4E265608D7F2

namespace mandos
{
/**
 * @brief The SimulationObject represents a physical object that can be simulated using Mandos.
 * It is composed of its degrees of freedom, represented by the physical state and its energies
 *
 * @tparam PhysicalStateT
 * Specifies the type of PhysicalState the SimulationObject has, like for example RigidBody dofs
 */
template <typename PhysicalStateT>
class SimulationObject
{
};

/**
 * @brief The energies supported by a SimulationObject<PhysicalStateT>

 * Some energies have sense only for a particular type of SimulationObject. For example, an object made of 3D particles
 * may have an elastic energy, like StVK. This same energy, doesn't make sense for a rigid body simulation object.
 *
 * Energies<> is specialized for each SimulationObject<> based on its PhysicalState.
 *
 * @tparam PhysicaStateT PhysicalState of the SimulationObject supporting this set of energies
 */
template <typename PhysicaStateT>
struct Energies {
};
}  // namespace mandos

#endif /* EA5CC9BC_9B26_421C_B4CB_4E265608D7F2 */
