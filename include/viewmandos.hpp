#ifndef VIEWMANDOS_H_
#define VIEWMANDOS_H_

#include "../src/linear_algebra.hpp"
#include "mandos.hpp"
#include "../src/mesh.hpp"
#include "../src/physics_state.hpp"
#include "../src/memory_pool.hpp"
#include "../src/simulation.hpp"
#include <vector>

enum KeyboardKeys {
    Key_NULL            = 0,        // Key: NULL, used for no key pressed
    // Alphanumeric keys
    Key_APOSTROPHE      = 39,       // Key: '
    Key_COMMA           = 44,       // Key: ,
    Key_MINUS           = 45,       // Key: -
    Key_PERIOD          = 46,       // Key: .
    Key_SLASH           = 47,       // Key: /
    Key_ZERO            = 48,       // Key: 0
    Key_ONE             = 49,       // Key: 1
    Key_TWO             = 50,       // Key: 2
    Key_THREE           = 51,       // Key: 3
    Key_FOUR            = 52,       // Key: 4
    Key_FIVE            = 53,       // Key: 5
    Key_SIX             = 54,       // Key: 6
    Key_SEVEN           = 55,       // Key: 7
    Key_EIGHT           = 56,       // Key: 8
    Key_NINE            = 57,       // Key: 9
    Key_SEMICOLON       = 59,       // Key: ;
    Key_EQUAL           = 61,       // Key: =
    Key_A               = 65,       // Key: A | a
    Key_B               = 66,       // Key: B | b
    Key_C               = 67,       // Key: C | c
    Key_D               = 68,       // Key: D | d
    Key_E               = 69,       // Key: E | e
    Key_F               = 70,       // Key: F | f
    Key_G               = 71,       // Key: G | g
    Key_H               = 72,       // Key: H | h
    Key_I               = 73,       // Key: I | i
    Key_J               = 74,       // Key: J | j
    Key_K               = 75,       // Key: K | k
    Key_L               = 76,       // Key: L | l
    Key_M               = 77,       // Key: M | m
    Key_N               = 78,       // Key: N | n
    Key_O               = 79,       // Key: O | o
    Key_P               = 80,       // Key: P | p
    Key_Q               = 81,       // Key: Q | q
    Key_R               = 82,       // Key: R | r
    Key_S               = 83,       // Key: S | s
    Key_T               = 84,       // Key: T | t
    Key_U               = 85,       // Key: U | u
    Key_V               = 86,       // Key: V | v
    Key_W               = 87,       // Key: W | w
    Key_X               = 88,       // Key: X | x
    Key_Y               = 89,       // Key: Y | y
    Key_Z               = 90,       // Key: Z | z
    Key_LEFT_BRACKET    = 91,       // Key: [
    Key_BACKSLASH       = 92,       // Key: '\'
    Key_RIGHT_BRACKET   = 93,       // Key: ]
    Key_GRAVE           = 96,       // Key: `
    // Function keys
    Key_SPACE           = 32,       // Key: Space
    Key_ESCAPE          = 256,      // Key: Esc
    Key_ENTER           = 257,      // Key: Enter
    Key_TAB             = 258,      // Key: Tab
    Key_BACKSPACE       = 259,      // Key: Backspace
    Key_INSERT          = 260,      // Key: Ins
    Key_DELETE          = 261,      // Key: Del
    Key_RIGHT           = 262,      // Key: Cursor right
    Key_LEFT            = 263,      // Key: Cursor left
    Key_DOWN            = 264,      // Key: Cursor down
    Key_UP              = 265,      // Key: Cursor up
    Key_PAGE_UP         = 266,      // Key: Page up
    Key_PAGE_DOWN       = 267,      // Key: Page down
    Key_HOME            = 268,      // Key: Home
    Key_END             = 269,      // Key: End
    Key_CAPS_LOCK       = 280,      // Key: Caps lock
    Key_SCROLL_LOCK     = 281,      // Key: Scroll down
    Key_NUM_LOCK        = 282,      // Key: Num lock
    Key_PRINT_SCREEN    = 283,      // Key: Print screen
    Key_PAUSE           = 284,      // Key: Pause
    Key_F1              = 290,      // Key: F1
    Key_F2              = 291,      // Key: F2
    Key_F3              = 292,      // Key: F3
    Key_F4              = 293,      // Key: F4
    Key_F5              = 294,      // Key: F5
    Key_F6              = 295,      // Key: F6
    Key_F7              = 296,      // Key: F7
    Key_F8              = 297,      // Key: F8
    Key_F9              = 298,      // Key: F9
    Key_F10             = 299,      // Key: F10
    Key_F11             = 300,      // Key: F11
    Key_F12             = 301,      // Key: F12
    Key_LEFT_SHIFT      = 340,      // Key: Shift left
    Key_LEFT_CONTROL    = 341,      // Key: Control left
    Key_LEFT_ALT        = 342,      // Key: Alt left
    Key_LEFT_SUPER      = 343,      // Key: Super left
    Key_RIGHT_SHIFT     = 344,      // Key: Shift right
    Key_RIGHT_CONTROL   = 345,      // Key: Control right
    Key_RIGHT_ALT       = 346,      // Key: Alt right
    Key_RIGHT_SUPER     = 347,      // Key: Super right
    Key_KB_MENU         = 348,      // Key: KB menu
    // Keypad keys
    Key_KP_0            = 320,      // Key: Keypad 0
    Key_KP_1            = 321,      // Key: Keypad 1
    Key_KP_2            = 322,      // Key: Keypad 2
    Key_KP_3            = 323,      // Key: Keypad 3
    Key_KP_4            = 324,      // Key: Keypad 4
    Key_KP_5            = 325,      // Key: Keypad 5
    Key_KP_6            = 326,      // Key: Keypad 6
    Key_KP_7            = 327,      // Key: Keypad 7
    Key_KP_8            = 328,      // Key: Keypad 8
    Key_KP_9            = 329,      // Key: Keypad 9
    Key_KP_DECIMAL      = 330,      // Key: Keypad .
    Key_KP_DIVIDE       = 331,      // Key: Keypad /
    Key_KP_MULTIPLY     = 332,      // Key: Keypad *
    Key_KP_SUBTRACT     = 333,      // Key: Keypad -
    Key_KP_ADD          = 334,      // Key: Keypad +
    Key_KP_ENTER        = 335,      // Key: Keypad Enter
    Key_KP_EQUAL        = 336,      // Key: Keypad =
};

struct MeshGPU {
    MeshGPU(const RenderMesh& mesh);
    ~MeshGPU();
    void updateData(const RenderMesh& mesh);
    int nVertices;
    unsigned int verticesVBO, texcoordsVBO, normalsVBO, tangentsVBO, VAO;
    float *vertices, *texcoords, *normals, *tangents;
};

struct MandosViewer {
    const int initialScreenWidth = 1600;
    const int initialScreenHeight = 900;
    MemoryPool mem_pool = MemoryPool(static_cast<size_t>(1e6) * sizeof(float));

    MandosViewer();
    MandosViewer(const Simulation* simulation);
    ~MandosViewer();

    void disable_render_logs();

    bool window_should_close();

    void begin_drawing();

    void end_drawing();

    bool is_key_pressed(int Key);

    void draw_particle(const ParticleHandle& particle, const PhysicsState& state);

    void draw_rigid_body(const RigidBodyHandle& rb, const PhysicsState& state, const MeshGPU& mesh);

    void draw_particles(const Simulation& simulation, const PhysicsState& state);

    void draw_rigid_bodies(const Simulation& simulation, const PhysicsState& state);

    void draw_FEM(const FEMHandle& fem, const PhysicsState& state, MeshGPU& gpuMesh, RenderMesh& renderMesh, SimulationMesh& simMesh);

    void draw_MassSpring(const MassSpringHandle& mass_spring, const PhysicsState& state, MeshGPU& gpuMesh, RenderMesh& renderMesh, SimulationMesh& simMesh);

    void draw_mesh(const Mat4& transform, const MeshGPU& mesh);

    void draw_simulation_state(const Simulation& simulation, const PhysicsState& state);

    void draw_springs_lines(const Simulation& simulation, const PhysicsState& state);

    void draw_springs(const Simulation& simulation, const PhysicsState& state);

    void draw_FEM_tetrahedrons_lines(const Simulation& simulation, const PhysicsState& state);

    void draw_FEM_tetrahedrons(const Simulation& simulation, const PhysicsState& state);

    void draw_particle_indices(const Simulation& simulation, const PhysicsState& state);

    void draw_vector(const Vec3& vector, const Vec3& origin);

    void draw_rods(const Simulation& simulation, const PhysicsState& state);

    void draw_colliders(const Simulation& simulation);

private:
    void initialize_graphical_context();
    void drawGUI();
    void drawSimulationVisualizationWindow();

    bool enable_normal_smoothing = false;
    bool enable_transparent_background = false;
    bool enable_draw_particle_indices = false;
    bool enable_draw_simulable_meshes = true;
    bool enable_draw_particles = false;
    bool enable_draw_rigid_bodies = false;
    bool enable_draw_springs = false;
    bool enable_draw_colliders = true;
    enum FEM_TETRAHEDRON_VISUALIZATION {TET_NONE, TET_MESH, TET_LINES};
    int enable_draw_fem_tetrahedrons = TET_NONE;

    bool enable_slice_plane = false;
    Vec4 slicePlane = Vec4(0.0f, 1.0f, 0.0f, 0.0f);
    bool slicePlaneGuizmoToggle = false;

    float simulationViewerWidth = static_cast<float>(initialScreenWidth);
    float simulationViewerHeight = static_cast<float>(initialScreenHeight);

    const Simulation* SavedSim = nullptr;
    const PhysicsState* SavedState = nullptr;
};

#endif // VIEWMANDOS_H_
