cmake_minimum_required(VERSION 3.7)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Obtain a file compile_commands.json used by ccls (through the plugin coc.nvim)
# to provide code completion in Neovim
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")

set(CMAKE_CXX_STANDARD 17)

# Name of the project
project(project)

# Find the OpenCV package
find_package( OpenCV REQUIRED )

# Create a variable with all the include directories
set(INCLUDE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/src

    # Include directories for the linked libraries
    ${OpenCV_INCLUDE_DIRS}
   )

include_directories(${INCLUDE})

# Create a variable with a link to all cpp files to compile
set(SOURCES
    main.cpp
    ${PROJECT_SOURCE_DIR}/src/NavierStokesSolver.cpp
    ${PROJECT_SOURCE_DIR}/src/boundaryConditions.cpp
   )

add_executable(main ${SOURCES})

# Link to the external libraries
target_link_libraries( main ${OpenCV_LIBS} )

# Get rid of the cmake_install.cmake file created
set(CMAKE_SKIP_INSTALL_RULES True)

# Move the file compile_commands.json to the directory of the project.
execute_process(COMMAND cp compile_commands.json ..) 
