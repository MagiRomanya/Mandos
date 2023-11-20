#include <iostream>
#include <raylib.h>
#include <rlgl.h>
#include <imgui.h>
#include <raymath.h>

#include "raylib_imgui.hpp"
#include "editor/edit_mode.hpp"

int main(int argc, char *argv[]) {
    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 1600;
    const int screenHeight = 900;

    InitWindow(screenWidth, screenHeight, "Mandos GUI");
    ImGui_initialize();
    auto io = ImGui::GetIO();

    edit_mode_loop();

    // De-Initialization
    //--------------------------------------------------------------------------------------
    ImGuiDeinitialize();

    CloseWindow();        // Close glfwWindow and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}
