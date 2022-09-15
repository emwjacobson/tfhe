/*
 * Free FFT and convolution (C++)
 *
 * Copyright (c) 2021 Project Nayuki. (MIT License)
 * https://www.nayuki.io/page/free-small-fft-in-multiple-languages
 *
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

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string.h>
#include <math.h>
#include <algorithm>
#include "fpga_constants.h"

// Will do the actual FFT multiplication
// Templated through SIZE and WINDOW
// eg. SIZE=2048 WINDOW=256 will perform an FFT stage on an array of size 2048 in 256 sized windows.
//     Loop `fft_mult_window` will run 8 times. Loop `fft_mult_mult` will run (256/2)=128 times.
//		 In theory `fft_mult_window` could be parallelized by a factor of 8, as each iteration of `fft_mult_window`
//		 does not depend on another. (Though it will depend on the array_partition factor)
// eg. SIZE=2048 WINDOW=2048 will perform an FFT stage on an array of size 2048 in 2048 sized windows.
//		 Loop `fft_mult_window` will run 1 time. Loop `fft_mult_mult` will run (2048/2)=1024 times.
//		 This example cannot be parallelized as the entire array is needed for the calculation.
template <size_t SIZE, size_t WINDOW>
void fft_mult(double *real_out, double *imag_out, const double *real_in, const double *imag_in) {
	constexpr size_t halfsize = WINDOW / 2;
	constexpr size_t tablestep = param_2N / WINDOW;

	fft_mult_window: for (size_t i = 0; i < SIZE; i += WINDOW) {
		fft_mult_mult: for (size_t j = 0; j < halfsize; j++) {
			#pragma HLS pipeline II=1
			const size_t first = j + i;
			const size_t second = j + i + halfsize;
			const int k = j * tablestep;
			// Load
			double tprel = real_in[second];
			double tpiml = imag_in[second];
			double tprej = real_in[first];
			double tpimj = imag_in[first];
			double tpcos = cosTable[k];
			double tpsin = sinTable[k];
			// Calc
			double calcre =  tprel * tpcos + tpiml * tpsin;
			double calcim = -tprel * tpsin + tpiml * tpcos;
			// Store
			real_out[second] = tprej - calcre;
			imag_out[second] = tpimj - calcim;
			real_out[first] = tprej + calcre;
			imag_out[first] = tpimj + calcim;
		}
	}
}

// Do FFT on SIZE=2048 with WINDOW=1024
// Can parallelize it as 2 SIZE=1024, WINDOW=1024
template <>
void fft_mult<2048, 1024>(double *real_out, double *imag_out, const double *real_in, const double *imag_in) {
	#pragma HLS dataflow
	fft_mult<1024, 1024>(&real_out[1024*0], &imag_out[1024*0], &real_in[1024*0], &imag_in[1024*0]);
	fft_mult<1024, 1024>(&real_out[1024*1], &imag_out[1024*1], &real_in[1024*1], &imag_in[1024*1]);
}

// Do FFT on SIZE=2048 with WINDOW=512
// Can parallelize it 4 ways with SIZE=512, WINDOW=512
template <>
void fft_mult<2048, 512>(double *real_out, double *imag_out, const double *real_in, const double *imag_in) {
	#pragma HLS dataflow
	fft_mult<512, 512>(&real_out[512*0], &imag_out[512*0], &real_in[512*0], &imag_in[512*0]);
	fft_mult<512, 512>(&real_out[512*1], &imag_out[512*1], &real_in[512*1], &imag_in[512*1]);
	fft_mult<512, 512>(&real_out[512*2], &imag_out[512*2], &real_in[512*2], &imag_in[512*2]);
	fft_mult<512, 512>(&real_out[512*3], &imag_out[512*3], &real_in[512*3], &imag_in[512*3]);
}

// Do FFT on SIZE=2048 with WINDOW=256
// Can parallelize it 8 ways with SIZE=256, WINDOW=256
template <>
void fft_mult<2048, 256>(double *real_out, double *imag_out, const double *real_in, const double *imag_in) {
	#pragma HLS dataflow
	fft_mult<256, 256>(&real_out[256*0], &imag_out[256*0], &real_in[256*0], &imag_in[256*0]);
	fft_mult<256, 256>(&real_out[256*1], &imag_out[256*1], &real_in[256*1], &imag_in[256*1]);
	fft_mult<256, 256>(&real_out[256*2], &imag_out[256*2], &real_in[256*2], &imag_in[256*2]);
	fft_mult<256, 256>(&real_out[256*3], &imag_out[256*3], &real_in[256*3], &imag_in[256*3]);
	fft_mult<256, 256>(&real_out[256*4], &imag_out[256*4], &real_in[256*4], &imag_in[256*4]);
	fft_mult<256, 256>(&real_out[256*5], &imag_out[256*5], &real_in[256*5], &imag_in[256*5]);
	fft_mult<256, 256>(&real_out[256*6], &imag_out[256*6], &real_in[256*6], &imag_in[256*6]);
	fft_mult<256, 256>(&real_out[256*7], &imag_out[256*7], &real_in[256*7], &imag_in[256*7]);
}

// Do FFT on SIZE=2048 with WINDOW=128
// Can parallelize it 8 ways (limited by array_partition) with SIZE=256, WINDOW=128
// THIS DOES NOT PASS DATAFLOW
template <>
void fft_mult<2048, 128>(double *real_out, double *imag_out, const double *real_in, const double *imag_in) {
	#pragma HLS dataflow
	fft_mult<256, 128>(&real_out[256*0], &imag_out[256*0], &real_in[256*0], &imag_in[256*0]);
	fft_mult<256, 128>(&real_out[256*1], &imag_out[256*1], &real_in[256*1], &imag_in[256*1]);
	fft_mult<256, 128>(&real_out[256*2], &imag_out[256*2], &real_in[256*2], &imag_in[256*2]);
	fft_mult<256, 128>(&real_out[256*3], &imag_out[256*3], &real_in[256*3], &imag_in[256*3]);
	fft_mult<256, 128>(&real_out[256*4], &imag_out[256*4], &real_in[256*4], &imag_in[256*4]);
	fft_mult<256, 128>(&real_out[256*5], &imag_out[256*5], &real_in[256*5], &imag_in[256*5]);
	fft_mult<256, 128>(&real_out[256*6], &imag_out[256*6], &real_in[256*6], &imag_in[256*6]);
	fft_mult<256, 128>(&real_out[256*7], &imag_out[256*7], &real_in[256*7], &imag_in[256*7]);
}

extern "C" {
	void fft_transform(double *real, double *imag) {
		double re1[param_2N];
		double im1[param_2N];
		double re2[param_2N];
		double im2[param_2N];

		#pragma HLS array_partition variable=re1 block factor=8 dim=1
		#pragma HLS array_partition variable=im1 block factor=8 dim=1
		#pragma HLS array_partition variable=re2 block factor=8 dim=1
		#pragma HLS array_partition variable=im2 block factor=8 dim=1

		// Bit reversal
		fft_bit_reverse: for (int i = 0; i < param_2N; i++) {
			uint64_t j = bit_reversed[i];
			re1[i] = real[j];
			im1[i] = imag[j];
		}

		// At this point, temp arrays hold the bit-reversed elements
		// Now we swap between the 2 memories to perform radix2fft
		// To do full FFT we need to do fft_mult for all powers of 2 from 2 to 2N.
		fft_mult<2048, 2>(re2, im2, re1, im1); // (re2,im2) = radix2fft(re1, im1)
		fft_mult<2048, 4>(re1, im1, re2, im2); // (re1,im1) = radix2fft(re2, im2)
		fft_mult<2048, 8>(re2, im2, re1, im1); // ...
		fft_mult<2048, 16>(re1, im1, re2, im2);
		fft_mult<2048, 32>(re2, im2, re1, im1);
		fft_mult<2048, 64>(re1, im1, re2, im2);
		fft_mult<2048, 128>(re2, im2, re1, im1);
		fft_mult<2048, 256>(re1, im1, re2, im2);
		fft_mult<2048, 512>(re2, im2, re1, im1);
		fft_mult<2048, 1024>(re1, im1, re2, im2);
		fft_mult<2048, 2048>(real, imag, re1, im1); // Final result stored in real, imag
	}
}