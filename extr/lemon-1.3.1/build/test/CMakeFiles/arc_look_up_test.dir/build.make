# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/antala/prog/lemon-1.3.1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/antala/prog/lemon-1.3.1/build

# Include any dependencies generated for this target.
include test/CMakeFiles/arc_look_up_test.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/arc_look_up_test.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/arc_look_up_test.dir/flags.make

test/CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.o: test/CMakeFiles/arc_look_up_test.dir/flags.make
test/CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.o: ../test/arc_look_up_test.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/antala/prog/lemon-1.3.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.o"
	cd /home/antala/prog/lemon-1.3.1/build/test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.o -c /home/antala/prog/lemon-1.3.1/test/arc_look_up_test.cc

test/CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.i"
	cd /home/antala/prog/lemon-1.3.1/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/antala/prog/lemon-1.3.1/test/arc_look_up_test.cc > CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.i

test/CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.s"
	cd /home/antala/prog/lemon-1.3.1/build/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/antala/prog/lemon-1.3.1/test/arc_look_up_test.cc -o CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.s

test/CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.o.requires:

.PHONY : test/CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.o.requires

test/CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.o.provides: test/CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.o.requires
	$(MAKE) -f test/CMakeFiles/arc_look_up_test.dir/build.make test/CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.o.provides.build
.PHONY : test/CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.o.provides

test/CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.o.provides.build: test/CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.o


# Object files for target arc_look_up_test
arc_look_up_test_OBJECTS = \
"CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.o"

# External object files for target arc_look_up_test
arc_look_up_test_EXTERNAL_OBJECTS =

test/arc_look_up_test: test/CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.o
test/arc_look_up_test: test/CMakeFiles/arc_look_up_test.dir/build.make
test/arc_look_up_test: lemon/libemon.a
test/arc_look_up_test: /usr/local/lib/libglpk.so
test/arc_look_up_test: test/CMakeFiles/arc_look_up_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/antala/prog/lemon-1.3.1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable arc_look_up_test"
	cd /home/antala/prog/lemon-1.3.1/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/arc_look_up_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/arc_look_up_test.dir/build: test/arc_look_up_test

.PHONY : test/CMakeFiles/arc_look_up_test.dir/build

test/CMakeFiles/arc_look_up_test.dir/requires: test/CMakeFiles/arc_look_up_test.dir/arc_look_up_test.cc.o.requires

.PHONY : test/CMakeFiles/arc_look_up_test.dir/requires

test/CMakeFiles/arc_look_up_test.dir/clean:
	cd /home/antala/prog/lemon-1.3.1/build/test && $(CMAKE_COMMAND) -P CMakeFiles/arc_look_up_test.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/arc_look_up_test.dir/clean

test/CMakeFiles/arc_look_up_test.dir/depend:
	cd /home/antala/prog/lemon-1.3.1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/antala/prog/lemon-1.3.1 /home/antala/prog/lemon-1.3.1/test /home/antala/prog/lemon-1.3.1/build /home/antala/prog/lemon-1.3.1/build/test /home/antala/prog/lemon-1.3.1/build/test/CMakeFiles/arc_look_up_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/arc_look_up_test.dir/depend

