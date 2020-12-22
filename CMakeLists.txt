cmake_minimum_required(VERSION 3.0)

set(pluginName XBraidForUG4) #[[todo]]
project(UG_PLUGIN_${pluginName})

set(UG4_XBraid_plugin_VERSION_MAJOR 0.1)
set(UG4_XBraid_plugin_VERSION_MINOR 0.0)

#[[ ----------  ----------  ---------- ---------- Compiling libraries ---------- ----------  ----------  ---------- ]]
if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/libs/xbraid/braid/libbraid.a")
    add_custom_command(OUTPUT braidlib
            COMMAND make braid
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/libs/xbraid
            )
else ()
     add_custom_command(
            OUTPUT braidlib
            COMMAND touch hello.txt
        )
endif ()
add_custom_target(run ALL DEPENDS braidlib)


#[[
IF(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -3")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
]]
set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS_INIT} -fPIC")



set(SOURCES
        src/xbraidforug4_plugin.cpp
        src/tools/Telvanni.h
        src/tools/Talasma.h
        src/PBraid.h
        src/BraidVectorStruct.h


        src/ITSGFBraidApp.h
        src/UG4BraidApp.h
        src/RGFBraidApp.h
        src/MGFBraidApp.h



        src/PBraidApp.h
        src/SpaceTimeCommunicator.h
        src/Scriptor.h


	    src/MemoryObserver.h
        )


if(UG_ROOT_CMAKE_PATH)
    MESSAGE("UG Root was found")
else(UG_ROOT_CMAKE_PATH)
    set(UG_ROOT_CMAKE_PATH /opt/ug/ugcore/cmake)
endif(UG_ROOT_CMAKE_PATH)

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




