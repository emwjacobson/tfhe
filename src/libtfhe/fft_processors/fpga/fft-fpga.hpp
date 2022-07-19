#ifndef __FFT_FPGA_HPP__
#define __FFT_FPGA_HPP__

#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY 1
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <CL/cl2.hpp>

std::vector<cl::Device> get_xilinx_devices();
char* read_binary_file(const std::string &xclbin_file_name, unsigned &nb);

struct FftTables {
	uint64_t n;
	uint64_t *bit_reversed;
	double *trig_tables;
};

class FPGACompute {
public:
  const int32_t _2N;
  const int32_t N;
  const int32_t Ns2;

  FPGACompute(const int32_t N);
  void fft_transform_reverse(const void *tables, double *real, double *imag);
private:
  cl::Context context;
  cl::CommandQueue q;
  cl::Kernel krnl_fft;

  void __fft_reverse_cpu(const void *tables, double *real, double *imag);
  void __fft_reverse_fpga(const void *tables, double *real, double *imag);
};

#endif