/*
 * Fast Fourier transform
 *
 * Copyright (c) 2016 Project Nayuki
 * https://www.nayuki.io/page/fast-fourier-transform-in-x86-assembly
 *
 * (MIT License)
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
 * the Software, and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 * - The Software is provided "as is", without warranty of any kind, express or
 *   implied, including but not limited to the warranties of merchantability,
 *   fitness for a particular purpose and noninfringement. In no event shall the
 *   authors or copyright holders be liable for any claim, damages or other
 *   liability, whether in an action of contract, tort or otherwise, arising from,
 *   out of or in connection with the Software or the use or other dealings in the
 *   Software.
 */

#ifndef FFT_FPGA_H
#define FFT_FPGA_H

#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY 1
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <CL/cl2.hpp>
#include <iostream>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>

#ifdef __cplusplus
extern "C" {
#endif

struct FftTables {
	uint64_t n;
	uint64_t *bit_reversed;
	double *trig_tables;
};

struct FftTables_Container {
  uint64_t n;
  uint64_t bit_reversed[2*1024 * sizeof(size_t)];
  double trig_tables[((2*1024) - 4) * 2 * sizeof(double)];
};

void* fft_init(size_t n);

void* fft_init_reverse(size_t n);

void fft_transform(const void *tables, double *real, double *imag);

void fft_transform_reverse(const void *tables, double *real, double *imag);

void fft_destroy(void *tables);

std::vector<cl::Device> get_xilinx_devices();

char* read_binary_file(const std::string &xclbin_file_name, unsigned &nb);

#ifdef __cplusplus
}
#endif

#endif //FFT_FPGA_H
