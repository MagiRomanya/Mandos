import os

from conan.tools.cmake import CMakeDeps, CMakeToolchain, cmake_layout
from conan.tools.env import VirtualRunEnv

from conan import ConanFile

# Force Conan 2.0
required_conan_version = ">=2.0.0"


class Mandos(ConanFile):
    name = "Mandos"
    settings = "os", "compiler", "build_type", "arch"

    def requirements(self):
        self.requires("eigen/3.4.0") # Algebra library
        self.requires("entt/3.7.1") # Used to provide ECS functionality
        self.requires("spdlog/1.14.1") # Logging
        self.requires("tracy/0.10") # Tracing performance
        self.requires("catch2/3.6.0") # Unit testing framework
    
    def configure(self):
        # No options to configure yet
        pass

    def generate(self):
        cmake = CMakeDeps(self)
        cmake.generate()

        toolchain = CMakeToolchain(self)
        toolchain.user_presets_path = "ConanPresets.json"
        toolchain.generate()

        venv = VirtualRunEnv(self)
        venv.generate()

    def layout(self):
        cmake_layout(self)
