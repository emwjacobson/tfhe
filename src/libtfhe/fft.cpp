#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "fft.h"
#include "tfhe_core.h"
#include "fpga.h"
#include "polynomials.h"

// M_PI isn't defined with C99 for whatever reason
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

/*---- Function implementations ----*/

// This is a C implementation that models the x86-64 AVX implementation.
void fft_transform(double *real, double *imag) {
	const uint64_t n = N_Values._2N;

	cl::Buffer real_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(double) * n);
	cl::Buffer imag_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(double) * n);

	double *real_map = (double *)fpga.q.enqueueMapBuffer(real_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(double) * n);
	double *imag_map = (double *)fpga.q.enqueueMapBuffer(imag_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(double) * n);

	memcpy(real_map, real, sizeof(double) * n);
	memcpy(imag_map, imag, sizeof(double) * n);

	fpga.k_fft_transform.setArg(0, real_buf);
	fpga.k_fft_transform.setArg(1, imag_buf);

	fpga.q.enqueueMigrateMemObjects({ real_buf, imag_buf }, 0 /* 0 means from host*/);
	fpga.q.enqueueTask(fpga.k_fft_transform);
	fpga.q.enqueueMigrateMemObjects({ real_buf, imag_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	fpga.q.finish();

	memcpy(real, real_map, sizeof(double) * n);
	memcpy(imag, imag_map, sizeof(double) * n);
}

// This is a C implementation that models the x86-64 AVX implementation.
void fft_transform_reverse(double *real, double *imag) {
	const uint64_t n = N_Values._2N;

	cl::Buffer real_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(double) * n);
	cl::Buffer imag_buf(fpga.context, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_WRITE, sizeof(double) * n);

	double *real_map = (double *)fpga.q.enqueueMapBuffer(real_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(double) * n);
	double *imag_map = (double *)fpga.q.enqueueMapBuffer(imag_buf, CL_TRUE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(double) * n);

	memcpy(real_map, real, sizeof(double) * n);
	memcpy(imag_map, imag, sizeof(double) * n);

	fpga.k_fft_transform_reverse.setArg(0, real_buf);
	fpga.k_fft_transform_reverse.setArg(1, imag_buf);

	fpga.q.enqueueMigrateMemObjects({ real_buf, imag_buf }, 0 /* 0 means from host*/);
	fpga.q.enqueueTask(fpga.k_fft_transform_reverse);
	fpga.q.enqueueMigrateMemObjects({ real_buf, imag_buf }, CL_MIGRATE_MEM_OBJECT_HOST);

	fpga.q.finish();

	memcpy(real, real_map, sizeof(double) * n);
	memcpy(imag, imag_map, sizeof(double) * n);
}
