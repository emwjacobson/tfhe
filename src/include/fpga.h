#ifndef __FPGA_H__
#define __FPGA_H__

#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY 1
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <CL/cl2.hpp>
#include "tfhe_core.h"

class FPGA_Processor {
public:
  FPGA_Processor();
  ~FPGA_Processor();
  cl::Context context;
  cl::CommandQueue q;
  cl::Kernel k_tfhe_blindRotate_FFT;

  cplx* omegaxminus1;

  cl::Buffer accum_buf;
  cl::Buffer bkFFT_buf;
  cl::Buffer bara_buf;

private:
  std::vector<cl::Device> get_xilinx_devices();
  std::vector<unsigned char> read_binary_file(const std::string &xclbin_file_name);

  cl::Program program;
};

extern FPGA_Processor fpga;

#endif