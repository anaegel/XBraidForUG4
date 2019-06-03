cmake_minimum_required(VERSION 3.0)

set(pluginName XBraidForUG4) #[[todo]]
project(UG_PLUGIN_${pluginName})

set(UG4_XBraid_plugin_VERSION_MAJOR 0.1)
set(UG4_XBraid_plugin_VERSION_MINOR 0.0)

#[[ ----------  ----------  ---------- ---------- Compiling libraries ---------- ----------  ----------  ---------- ]]

#[[ XBraid libraries TODO compile if not exists ]]
option(RECOMPILE_UPDATE "Update XBraid" On)
option(XBRAID_RECOMPILE "Recompile XBraid" On)


if (XBRAID_RECOMPILE OR NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/libs/xbraid/braid/libbraid.a")
    add_custom_command(OUTPUT hello.txt
            COMMAND make braid
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/libs/xbraid
            )
else ()
    #[[ add_custom_command(
            OUTPUT hello.txt
            COMMAND touch hello.txt
        ]]
endif ()
add_custom_target(run ALL DEPENDS hello.txt)

#[[
IF(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -3")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
]]
set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS_INIT} -fPIC")

set(SOURCES
        src/xbraidforug4_plugin.cpp
        src/stdafx.h
        src/BraidVectorStruct.h

        src/PBraidApp.h
        src/GFBraidApp.h
        src/RGFBraidApp.h
        src/PBraid.h

        src/XCommunicator.h
        )


if (TARGET run)
    set(UG_ROOT_CMAKE_PATH "/home/parnet/Projects/ug4/ugcore/cmake") #[[todo]]
endif ()
include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

if (buildEmbeddedPlugins)
    EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
else (buildEmbeddedPlugins)
    add_library(${pluginName} SHARED ${SOURCES})
    #add_library(${pluginName} SHARED ${CMAKE_CURRENT_SOURCE_DIR}/libs/xbraid/braid/libbraid.a)
    target_link_libraries(${pluginName} ${CMAKE_CURRENT_SOURCE_DIR}/libs/xbraid/braid/libbraid.a)
    target_link_libraries(${pluginName} ug4)
    #add_dependencies(${pluginName} UG_PLUGIN_ConvectionDiffusion)
    LINK_DIRECTORIES(${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    MESSAGE(${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    target_link_libraries(${pluginName} ConvectionDiffusion)
endif (buildEmbeddedPlugins)


#[[if(TARGET run)
    #[[ MPI Libraries
    find_package(MPI)
    include_directories(${MPI_INCLUDE_PATH})
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

    add_executable(plugin_xbraid ${SOURCES})
    add_dependencies(plugin_xbraid run)

    target_link_libraries(plugin_xbraid ${MPI_LIBRARIES})
    target_link_libraries(plugin_xbraid ${CMAKE_CURRENT_SOURCE_DIR}/libs/xbraid/braid/libbraid.a)
endif()]]




