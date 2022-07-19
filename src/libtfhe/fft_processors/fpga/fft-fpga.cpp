#include <vector>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <chrono>
#include <string>
#include "fft.h"
#include "fft-fpga.hpp"
#include <CL/cl2.hpp>


FPGACompute::FPGACompute(const int32_t N) : _2N(2*N),N(N),Ns2(N/2)
{
  cl_int err;
  unsigned fileBufSize;
  std::vector<cl::Device> devices = get_xilinx_devices();
  devices.resize(1);
  cl::Device device = devices[0];
  this->context = cl::Context(device, NULL, NULL, NULL, &err);
  char *fileBuf = read_binary_file("fft_fpga.xclbin", fileBufSize);
  cl::Program::Binaries bins{{fileBuf, fileBufSize}};
  cl::Program program(context, devices, bins, NULL, &err);
  this->q = cl::CommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &err);
  this->krnl_fft = cl::Kernel(program, "fft_transform_reverse", &err);
}

void FPGACompute::fft_transform_reverse(const void *tables, double *real, double *imag) {
}

void FPGACompute::__fft_reverse_cpu(const void *tables, double *real, double *imag) {

	struct FftTables *tbl = (struct FftTables *)tables;
	uint64_t n = tbl->n;

	// Bit-reversed addressing permutation
	uint64_t i;
	uint64_t *bitreversed = tbl->bit_reversed;
	for (i = 0; i < n; i++) {
		uint64_t j = bitreversed[i];
		if (i < j) {
			double tp0re = real[i];
			double tp0im = imag[i];
			double tp1re = real[j];
			double tp1im = imag[j];
			real[i] = tp1re;
			imag[i] = tp1im;
			real[j] = tp0re;
			imag[j] = tp0im;
		}
	}

	// Size 2 merge (special)
	if (n >= 2) {
		for (i = 0; i < n; i += 2) {
			double tpre = real[i];
			double tpim = imag[i];
			real[i] += real[i + 1];
			imag[i] += imag[i + 1];
			real[i + 1] = tpre - real[i + 1];
			imag[i + 1] = tpim - imag[i + 1];
		}
	}

	// Size 4 merge (special)
	if (n >= 4) {
		for (i = 0; i < n; i += 4) {
			// Even indices
			double tpre, tpim;
			tpre = real[i];
			tpim = imag[i];
			real[i] += real[i + 2];
			imag[i] += imag[i + 2];
			real[i + 2] = tpre - real[i + 2];
			imag[i + 2] = tpim - imag[i + 2];
			// Odd indices
			tpre = real[i + 1];
			tpim = imag[i + 1];
			real[i + 1] -= imag[i + 3];
			imag[i + 1] += real[i + 3];
			tpre += imag[i + 3];
			tpim -= real[i + 3];
			real[i + 3] = tpre;
			imag[i + 3] = tpim;
		}
	}

	// Size 8 and larger merges (general)
	double *trigtables = tbl->trig_tables;
	uint64_t size;
	for (size = 8; size <= n; size <<= 1) {
		uint64_t halfsize = size >> 1;
		uint64_t i;
		for (i = 0; i < n; i += size) {
			uint64_t j, off;
			for (j = 0, off = 0; j < halfsize; j += 4, off += 8) {
				uint64_t k;
				for (k = 0; k < 4; k++) {  // To simulate x86 AVX 4-vectors
					uint64_t vi = i + j + k;  // Vector index
					uint64_t ti = off + k;    // Table index
					double re = real[vi + halfsize];
					double im = imag[vi + halfsize];
					double tpre = re * trigtables[ti] + im * trigtables[ti + 4];
					double tpim = im * trigtables[ti] - re * trigtables[ti + 4];
					real[vi + halfsize] = real[vi] - tpre;
					imag[vi + halfsize] = imag[vi] - tpim;
					real[vi] += tpre;
					imag[vi] += tpim;
				}
			}
		}
		if (size == n)
			break;
		trigtables += size;
	}
}

void FPGACompute::__fft_reverse_fpga(const void *tables, double *real, double *imag) {
  struct FftTables *tbl = (struct FftTables *)tables;
	uint64_t n = tbl->n;

  // Size from `fft_init_reverse` in `fft-x8664-avx-aux.c`
  cl::Buffer tables_buf(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY, sizeof(FftTables) + (_2N * sizeof(size_t)) + ((_2N - 4) * 2 * sizeof(double)));
  cl::Buffer real_buf(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY, sizeof(double) * _2N);
  cl::Buffer imag_buf(context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY, sizeof(double) * _2N);

  FftTables *tables_mapped = (FftTables *)q.enqueueMapBuffer(tables_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(FftTables) + (_2N * sizeof(size_t)) + ((_2N - 4) * 2 * sizeof(double)));
  double *real_mapped = (double *)q.enqueueMapBuffer(tables_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(double) * _2N);
  double *imag_mapped = (double *)q.enqueueMapBuffer(tables_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(double) * _2N);

  tables_mapped->n = tbl->n;
  memcpy(tables_mapped->bit_reversed, tbl->bit_reversed, _2N * sizeof(size_t));
  memcpy(tables_mapped->trig_tables, tbl->trig_tables, (n - 4) * 2 * sizeof(double));
  for(int i = 0; i<_2N; i++) {
    real_mapped[i] = real[i];
    imag_mapped[i] = imag[i];
  }

  this->krnl_fft.setArg(0, tables_buf);
  this->krnl_fft.setArg(1, real_buf);
  this->krnl_fft.setArg(2, imag_buf);

  this->q.enqueueMigrateMemObjects({tables_buf, real_buf, imag_buf}, 0);
  this->q.enqueueTask(this->krnl_fft);
  this->q.enqueueMigrateMemObjects({tables_buf, real_buf, imag_buf}, CL_MIGRATE_MEM_OBJECT_HOST);

  tbl->n = tables_mapped->n;
  memcpy(tbl->bit_reversed, tables_mapped->bit_reversed, _2N * sizeof(size_t));
  memcpy(tbl->trig_tables, tables_mapped->trig_tables, (n - 4) * 2 * sizeof(double));
  for(int i = 0; i<_2N; i++) {
    real[i] = real_mapped[i];
    imag[i] = imag_mapped[i];
  }

  this->q.finish();
}












// Helper Functions

std::vector<cl::Device> get_xilinx_devices()
{
  size_t i;
  cl_int err;
  std::vector<cl::Platform> platforms;
  err = cl::Platform::get(&platforms);
  cl::Platform platform;
  for (i = 0; i < platforms.size(); i++)
  {
    platform = platforms[i];
    std::string platformName = platform.getInfo<CL_PLATFORM_NAME>(&err);
    if (platformName == "Xilinx")
    {
      std::cout << "INFO: Found Xilinx Platform" << std::endl;
      break;
    }
  }
  if (i == platforms.size())
  {
    std::cout << "ERROR: Failed to find Xilinx platform" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Getting ACCELERATOR Devices and selecting 1st such device
  std::vector<cl::Device> devices;
  err = platform.getDevices(CL_DEVICE_TYPE_ACCELERATOR, &devices);
  return devices;
}

char *read_binary_file(const std::string &xclbin_file_name, unsigned &nb)
{
  if (access(xclbin_file_name.c_str(), R_OK) != 0)
  {
    printf("ERROR: %s xclbin not available please build\n", xclbin_file_name.c_str());
    exit(EXIT_FAILURE);
  }
  // Loading XCL Bin into char buffer
  std::cout << "INFO: Loading '" << xclbin_file_name << "'\n";
  std::ifstream bin_file(xclbin_file_name.c_str(), std::ifstream::binary);
  bin_file.seekg(0, bin_file.end);
  nb = bin_file.tellg();
  bin_file.seekg(0, bin_file.beg);
  char *buf = new char[nb];
  bin_file.read(buf, nb);
  return buf;
}