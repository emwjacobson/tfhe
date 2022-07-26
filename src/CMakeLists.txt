cmake_minimum_required(VERSION 3.0)

#set(CMAKE_VERBOSE_MAKEFILE ON)
#set(CMAKE_RULE_MESSAGES ON)

set(CMAKE_CONFIGURATION_TYPES Debug Optim CACHE TYPE INTERNAL FORCE)
set(CMAKE_BUILD_TYPE "optim" CACHE STRING "Build Type: Debug or Optim")
set(ENABLE_FPGA ON CACHE BOOL "Enable the FPGA FFT processor")
set(ENABLE_TESTS OFF CACHE BOOL "Build the tests (requires googletest)")

project(tfhe)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "AppleClang")
  # https://stackoverflow.com/a/16229679
  set(CLANG_FLAGS "-stdlib=libc++")
endif()

# -std=c99 seems to be required in Travis tests for whatever reason
set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -std=c11")

set(CMAKE_CXX_FLAGS_DEBUG "${CLANG_FLAGS} -std=c++11 -g3 -O0 -Wall -Werror -lpthread -lrt -lstdc++")
set(CMAKE_C_FLAGS_DEBUG "-g3 -O0 -Wall -Werror -lpthread -lrt -lstdc++")

set(CMAKE_CXX_FLAGS_OPTIM "${CLANG_FLAGS} -std=c++11 -g3 -march=native -O2 -DNDEBUG -funroll-loops -Wall -Werror -lpthread -lrt -lstdc++")
set(CMAKE_C_FLAGS_OPTIM "-g3 -march=native -O3 -DNDEBUG -funroll-loops -Wall -Werror -lpthread -lrt -lstdc++")

set(CMAKE_CXX_FLAGS_RELEASE "${CLANG_FLAGS} -std=c++11 -g0 -march=native -O2 -DNDEBUG -funroll-loops -Wall -Werror -lpthread -lrt -lstdc++")
set(CMAKE_C_FLAGS_RELEASE "-g0 -march=native -O3 -DNDEBUG -funroll-loops -Wall -Werror -lpthread -lrt -lstdc++")

find_package(OpenCL REQUIRED)

include_directories("include")
file(GLOB TFHE_HEADERS include/*.h)

install(FILES ${TFHE_HEADERS}
    DESTINATION include/tfhe
    PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)

# include the lib and the tests
add_subdirectory(libtfhe)
if (ENABLE_TESTS)
enable_testing()
add_subdirectory(test)
endif (ENABLE_TESTS)
