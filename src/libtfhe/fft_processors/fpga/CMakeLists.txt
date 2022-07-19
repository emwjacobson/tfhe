cmake_minimum_required(VERSION 3.0)

set(SRCS_PORTABLE
    fft-x8664-avx-aux.c
    fft-model-of-x8664-avx.c
    fft_processor_nayuki.cpp
    lagrangehalfc_impl.cpp
    fft-fpga.cpp
    )

set(HEADERS
    fft-fpga.hpp
    fft.h
    lagrangehalfc_impl.h
    )

if (ENABLE_FPGA_FFT)
    add_library(tfhe-fft-fpga OBJECT ${SRCS_PORTABLE} ${HEADERS})
    set_property(TARGET tfhe-fft-fpga PROPERTY POSITION_INDEPENDENT_CODE ON)
endif (ENABLE_FPGA_FFT)
