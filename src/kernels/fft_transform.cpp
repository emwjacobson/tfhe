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

template <size_t SIZE>
void fft_mult(double *real_out, double *imag_out, const double *real_in, const double *imag_in) {
	constexpr size_t halfsize = SIZE / 2;
	constexpr size_t tablestep = param_2N / SIZE;

	for (size_t i = 0; i < param_2N; i += SIZE) {
		for (size_t j = 0; j < halfsize; j++) {
			// #pragma HLS dependence variable=real inter false // Amin recommended intra, inter seems to produce fewer errors?
			// #pragma HLS dependence variable=imag inter false
			// #pragma HLS dependence variable=real intra false
			// #pragma HLS dependence variable=imag intra false
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

extern "C" {
	void fft_transform(double *_real, double *_imag) {
		double re1[param_2N];
		double im1[param_2N];

		double re2[param_2N];
		double im2[param_2N];

		memcpy(re1, _real, sizeof(double) * param_2N);
		memcpy(im1, _imag, sizeof(double) * param_2N);

		// Bit reversal
		for (int i = 0; i < param_2N; i++) {
			// #pragma HLS dependence variable=real intra false
			// #pragma HLS dependence variable=imag intra false
			#pragma HLS pipeline II=1
			uint64_t j = bit_reversed[i];
			if (j > i) {
				re1[i] = _real[j];
				im1[i] = _imag[j];

				re1[j] = _real[i];
				im1[j] = _imag[i];
				// double tmpre = re1[i];
				// re1[i] = re1[j];
				// re1[j] = tmpre;

				// double tmpim = im1[i];
				// im1[i] = im1[j];
				// im1[j] = tmpim;
			}
		}

		// At this point, temp arrays hold the bit-reversed elements
		// Now we swap between the 2 memories to perform radix2fft
		fft_mult<2>(re2, im2, re1, im1); // (re2,im2) = radix2fft(re1, im1)
		fft_mult<4>(re1, im1, re2, im2); // (re1,im1) = radix2fft(re2,im2)
		fft_mult<8>(re2, im2, re1, im1); // ...
		fft_mult<16>(re1, im1, re2, im2);
		fft_mult<32>(re2, im2, re1, im1);
		fft_mult<64>(re1, im1, re2, im2);
		fft_mult<128>(re2, im2, re1, im1);
		fft_mult<256>(re1, im1, re2, im2);
		fft_mult<512>(re2, im2, re1, im1);
		fft_mult<1024>(re1, im1, re2, im2);
		fft_mult<2048>(re2, im2, re1, im1); // Final result stored in re2, im2

		memcpy(_real, re2, sizeof(double) * param_2N);
		memcpy(_imag, im2, sizeof(double) * param_2N);
	}
}