#include <iostream>
#include <raylib.h>
#include <rlgl.h>
#include <rlImGui.h>
#include <imgui.h>
#include <raymath.h>

#include "gui/edit_mode.hpp"

int main(int argc, char *argv[]) {
    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 1600;
    const int screenHeight = 900;

    InitWindow(screenWidth, screenHeight, "Mandos GUI");
    rlImGuiSetup(false);

    edit_mode_loop();

    // De-Initialization
    //--------------------------------------------------------------------------------------
    rlImGuiShutdown();
    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}
