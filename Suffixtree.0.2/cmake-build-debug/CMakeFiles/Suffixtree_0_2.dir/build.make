# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2020.1.1\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2020.1.1\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\AmirM-Negar\CLionProjects\Suffixtree.0.2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\AmirM-Negar\CLionProjects\Suffixtree.0.2\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Suffixtree_0_2.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Suffixtree_0_2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Suffixtree_0_2.dir/flags.make

CMakeFiles/Suffixtree_0_2.dir/main.cpp.obj: CMakeFiles/Suffixtree_0_2.dir/flags.make
CMakeFiles/Suffixtree_0_2.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\AmirM-Negar\CLionProjects\Suffixtree.0.2\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Suffixtree_0_2.dir/main.cpp.obj"
	C:\MinGW\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\Suffixtree_0_2.dir\main.cpp.obj -c C:\Users\AmirM-Negar\CLionProjects\Suffixtree.0.2\main.cpp

CMakeFiles/Suffixtree_0_2.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Suffixtree_0_2.dir/main.cpp.i"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\AmirM-Negar\CLionProjects\Suffixtree.0.2\main.cpp > CMakeFiles\Suffixtree_0_2.dir\main.cpp.i

CMakeFiles/Suffixtree_0_2.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Suffixtree_0_2.dir/main.cpp.s"
	C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\AmirM-Negar\CLionProjects\Suffixtree.0.2\main.cpp -o CMakeFiles\Suffixtree_0_2.dir\main.cpp.s

# Object files for target Suffixtree_0_2
Suffixtree_0_2_OBJECTS = \
"CMakeFiles/Suffixtree_0_2.dir/main.cpp.obj"

# External object files for target Suffixtree_0_2
Suffixtree_0_2_EXTERNAL_OBJECTS =

Suffixtree_0_2.exe: CMakeFiles/Suffixtree_0_2.dir/main.cpp.obj
Suffixtree_0_2.exe: CMakeFiles/Suffixtree_0_2.dir/build.make
Suffixtree_0_2.exe: C:/MinGW/lib/libsfml-graphics-d.a
Suffixtree_0_2.exe: C:/MinGW/lib/libsfml-window-d.a
Suffixtree_0_2.exe: C:/MinGW/lib/libsfml-system-d.a
Suffixtree_0_2.exe: CMakeFiles/Suffixtree_0_2.dir/linklibs.rsp
Suffixtree_0_2.exe: CMakeFiles/Suffixtree_0_2.dir/objects1.rsp
Suffixtree_0_2.exe: CMakeFiles/Suffixtree_0_2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\AmirM-Negar\CLionProjects\Suffixtree.0.2\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Suffixtree_0_2.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\Suffixtree_0_2.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Suffixtree_0_2.dir/build: Suffixtree_0_2.exe

.PHONY : CMakeFiles/Suffixtree_0_2.dir/build

CMakeFiles/Suffixtree_0_2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\Suffixtree_0_2.dir\cmake_clean.cmake
.PHONY : CMakeFiles/Suffixtree_0_2.dir/clean

CMakeFiles/Suffixtree_0_2.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\AmirM-Negar\CLionProjects\Suffixtree.0.2 C:\Users\AmirM-Negar\CLionProjects\Suffixtree.0.2 C:\Users\AmirM-Negar\CLionProjects\Suffixtree.0.2\cmake-build-debug C:\Users\AmirM-Negar\CLionProjects\Suffixtree.0.2\cmake-build-debug C:\Users\AmirM-Negar\CLionProjects\Suffixtree.0.2\cmake-build-debug\CMakeFiles\Suffixtree_0_2.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Suffixtree_0_2.dir/depend

