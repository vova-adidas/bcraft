cmake_minimum_required (VERSION 3.10)

set(CMAKE_C_COMPILER clang-cl.exe)
set(CMAKE_CXX_COMPILER clang-cl.exe)

project(bcraft)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/external")
find_package(SDL2 REQUIRED)
find_package(xsimd REQUIRED)
find_package(glad REQUIRED)
find_package(glm REQUIRED)
find_package(Stb REQUIRED)
if(NOT MSVC)
	find_package(OpenMP REQUIRED)
endif()

add_executable (bcraft
	"src/concurrency/parallel_for.hpp"
	"src/concurrency/concurrent_queue.hpp"
	"src/concurrency/thread_pool.hpp"
	"src/concurrency/spinlock.hpp"
	"src/concurrency/spinlock.hpp"
	"src/rendering/software/render_target.hpp"
	"src/rendering/software/render_constants.hpp"
	"src/rendering/software/warp.hpp"
	"src/rendering/software/software_renderer.hpp"
	"src/rendering/meshing.hpp" 
	"src/rendering/mesh_buffer_pool.hpp" 
	"src/rendering/view.hpp" 
	"src/utility/poskey.hpp"
	"src/core/world.hpp"
	"src/core/chunk.hpp" 
	"src/core/chunk_generation.hpp" 
	"src/core/chunk_lighting.hpp" 
	"src/main.cpp"
	src/core/chunk_generation.cpp
	src/core/chunk_lighting.cpp
	src/core/actor.hpp
	src/core/world.cpp
	src/concurrency/thread_pool.cpp
	src/display/FramebufferRenderer.hpp
	src/display/FramebufferRenderer.cpp
	src/rendering/software/texture_rgba.hpp
)

set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")

set_property(TARGET bcraft PROPERTY CXX_STANDARD 23)

if(MSVC)
	target_compile_definitions(bcraft PRIVATE BCRAFT_STDEXECUTION)
	target_compile_options(bcraft PRIVATE /O2 /Ot /Oi /GT /GL /Qvec- /arch:AVX2)
else()
	target_compile_definitions(bcraft PRIVATE BCRAFT_OPENMP)
	target_compile_options(bcraft PRIVATE -O2 -mavx2 -fno-vectorize -fno-slp-vectorize -flto)
endif()


set_property(TARGET bcraft PROPERTY CXX_STANDARD 23)
set_property(TARGET bcraft PROPERTY CMAKE_C_EXTENSIONS OFF)
set_property(TARGET bcraft PROPERTY CMAKE_CXX_EXTENSIONS OFF)
target_link_libraries(bcraft PRIVATE SDL2::SDL2 SDL2::SDL2main)
target_link_libraries(bcraft PRIVATE xsimd)
target_link_libraries(bcraft PRIVATE glm::glm)
target_link_libraries(bcraft PRIVATE glad::glad)
target_include_directories(bcraft PRIVATE ${Stb_INCLUDE_DIR})

if(NOT MSVC)
	target_link_libraries(bcraft PRIVATE OpenMP::OpenMP_CXX)
endif ()

set(TEXTURE_ATLAS "${CMAKE_CURRENT_SOURCE_DIR}/texture.png")

add_custom_command(
		TARGET bcraft POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy
		"${TEXTURE_ATLAS}"
		"$<TARGET_FILE_DIR:bcraft>/texture.png"
)
