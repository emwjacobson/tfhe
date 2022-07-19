cmake_minimum_required(VERSION 3.0)

set(SRCS
    fft_processor_fpga.cpp
    lagrangehalfc_impl.cpp
    fft.c
    )

set(HEADERS
    fft.h
    lagrangehalfc_impl.h
    )

if (ENABLE_FPGA)
    add_library(tfhe-fft-fpga OBJECT ${SRCS} ${HEADERS})
    set_property(TARGET tfhe-fft-fpga PROPERTY POSITION_INDEPENDENT_CODE ON)
endif (ENABLE_FPGA)
