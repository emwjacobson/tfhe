cmake_minimum_required(VERSION 3.0)

add_subdirectory(googletest EXCLUDE_FROM_ALL)

include_directories(${GTEST_INCLUDE_DIRS})

set(GOOGLETEST_SOURCES
        arithmetic_test.cpp
        lwe_test.cpp
        polynomial_test.cpp
        tlwe_test.cpp
        tgsw_test.cpp
        tlwe_fft_test.cpp
        tgsw_fft_test.cpp
        lwekeyswitch_test.cpp
        bootstrapping_test.cpp
        bootstrapping_test_fft.cpp
        io_test.cpp
        lagrangehalfc_test.cpp
        boots_gates_test.cpp
        fakes/lagrangehalfc.h
        fakes/lwe.h
        fakes/lwe-bootstrapping-fft.h
        fakes/lwe-bootstrapping.h
        fakes/tgsw-fft.h
        fakes/tgsw.h
        fakes/tlwe.h
        fakes/tlwe-fft.h
        fakes/lwe-keyswitch.h
        )

set(CPP_ITESTS
        test-bootstrapping-fft
        test-decomp-tgsw
        test-lwe
        test-multiplication
        test-tlwe
        test-gate-bootstrapping
        test-addition-boot
        test-long-run
        )

# set(C_ITESTS
#         test-c-binding
#         )

set(RUNTIME_LIBS
        tfhe-core
        ${OpenCL_LIBRARIES}
        )

# #the unittests are compiled with the google test framework
# add_executable(unittests-fpga ${GOOGLETEST_SOURCES} ${TFHE_HEADERS})
# target_link_libraries(unittests-fpga ${RUNTIME_LIBS} gtest gtest_main -lpthread -lpthread -lrt -lstdc++)
# add_test(unittests-fpga unittests-fpga)

#the integration tests must be single source code, and are compiled as a standalone application
#we first compile the C++ tests
foreach (CPP_ITEST ${CPP_ITESTS})
add_executable(${CPP_ITEST}-fpga ${CPP_ITEST}.cpp ${TFHE_HEADERS})
target_link_libraries(${CPP_ITEST}-fpga ${RUNTIME_LIBS} -lpthread -lrt -lstdc++)
endforeach (CPP_ITEST)

#then the C tests
foreach (C_ITEST ${C_ITESTS})
add_executable(${C_ITEST}-fpga ${C_ITEST}.c ${TFHE_HEADERS})
target_link_libraries(${C_ITEST}-fpga ${RUNTIME_LIBS})
endforeach (C_ITEST)
