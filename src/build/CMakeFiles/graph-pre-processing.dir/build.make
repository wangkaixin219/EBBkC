# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/kaixin/EBBkC/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/kaixin/EBBkC/src/build

# Include any dependencies generated for this target.
include CMakeFiles/graph-pre-processing.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/graph-pre-processing.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/graph-pre-processing.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/graph-pre-processing.dir/flags.make

CMakeFiles/graph-pre-processing.dir/truss/util/graph/graph_preprocessing.cpp.o: CMakeFiles/graph-pre-processing.dir/flags.make
CMakeFiles/graph-pre-processing.dir/truss/util/graph/graph_preprocessing.cpp.o: /home/kaixin/EBBkC/src/truss/util/graph/graph_preprocessing.cpp
CMakeFiles/graph-pre-processing.dir/truss/util/graph/graph_preprocessing.cpp.o: CMakeFiles/graph-pre-processing.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/kaixin/EBBkC/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/graph-pre-processing.dir/truss/util/graph/graph_preprocessing.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/graph-pre-processing.dir/truss/util/graph/graph_preprocessing.cpp.o -MF CMakeFiles/graph-pre-processing.dir/truss/util/graph/graph_preprocessing.cpp.o.d -o CMakeFiles/graph-pre-processing.dir/truss/util/graph/graph_preprocessing.cpp.o -c /home/kaixin/EBBkC/src/truss/util/graph/graph_preprocessing.cpp

CMakeFiles/graph-pre-processing.dir/truss/util/graph/graph_preprocessing.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/graph-pre-processing.dir/truss/util/graph/graph_preprocessing.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/kaixin/EBBkC/src/truss/util/graph/graph_preprocessing.cpp > CMakeFiles/graph-pre-processing.dir/truss/util/graph/graph_preprocessing.cpp.i

CMakeFiles/graph-pre-processing.dir/truss/util/graph/graph_preprocessing.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/graph-pre-processing.dir/truss/util/graph/graph_preprocessing.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/kaixin/EBBkC/src/truss/util/graph/graph_preprocessing.cpp -o CMakeFiles/graph-pre-processing.dir/truss/util/graph/graph_preprocessing.cpp.s

# Object files for target graph-pre-processing
graph__pre__processing_OBJECTS = \
"CMakeFiles/graph-pre-processing.dir/truss/util/graph/graph_preprocessing.cpp.o"

# External object files for target graph-pre-processing
graph__pre__processing_EXTERNAL_OBJECTS =

libgraph-pre-processing.a: CMakeFiles/graph-pre-processing.dir/truss/util/graph/graph_preprocessing.cpp.o
libgraph-pre-processing.a: CMakeFiles/graph-pre-processing.dir/build.make
libgraph-pre-processing.a: CMakeFiles/graph-pre-processing.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/kaixin/EBBkC/src/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libgraph-pre-processing.a"
	$(CMAKE_COMMAND) -P CMakeFiles/graph-pre-processing.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/graph-pre-processing.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/graph-pre-processing.dir/build: libgraph-pre-processing.a
.PHONY : CMakeFiles/graph-pre-processing.dir/build

CMakeFiles/graph-pre-processing.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/graph-pre-processing.dir/cmake_clean.cmake
.PHONY : CMakeFiles/graph-pre-processing.dir/clean

CMakeFiles/graph-pre-processing.dir/depend:
	cd /home/kaixin/EBBkC/src/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/kaixin/EBBkC/src /home/kaixin/EBBkC/src /home/kaixin/EBBkC/src/build /home/kaixin/EBBkC/src/build /home/kaixin/EBBkC/src/build/CMakeFiles/graph-pre-processing.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/graph-pre-processing.dir/depend

