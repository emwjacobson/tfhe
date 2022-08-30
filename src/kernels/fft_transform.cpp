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

extern "C" {
	void fft_transform(double *real, double *imag) {
		// Bit reversal
		for (int i = 0; i < param_2N; i++) {
			uint64_t j = bit_reversed[i];
			if (j > i) {
				double tmpre = real[i];
				real[i] = real[j];
				real[j] = tmpre;

				double tmpim = imag[i];
				imag[i] = imag[j];
				imag[j] = tmpim;
			}
		}

		// size = 2 (1024 parallel)
		for (size_t i = 0; i < param_2N; i += 2) {
			#pragma HLS unroll factor=16
			int k = 0;
			for (size_t j = i; j < i + 1; j++) {
				#pragma HLS pipeline II=1
				size_t l = j + 1;
				double tpre =  real[l] * cosTable[k] + imag[l] * sinTable[k];
				double tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k];
				real[l] = real[j] - tpre;
				imag[l] = imag[j] - tpim;
				real[j] += tpre;
				imag[j] += tpim;
				k += 1024;
			}
		}

		// size = 4 (512 parallel)
		for (size_t i = 0; i < param_2N; i += 4) {
			#pragma HLS unroll factor=8
			int k = 0;
			for (size_t j = i; j < i + 2; j++) {
				#pragma HLS pipeline II=1
				size_t l = j + 2;
				double tpre =  real[l] * cosTable[k] + imag[l] * sinTable[k];
				double tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k];
				real[l] = real[j] - tpre;
				imag[l] = imag[j] - tpim;
				real[j] += tpre;
				imag[j] += tpim;
				k += 512;
			}
		}

		// size = 8 (256 parallel)
		for (size_t i = 0; i < param_2N; i += 8) {
			#pragma HLS unroll factor=4
			int k = 0;
			for (size_t j = i; j < i + 4; j++) {
				#pragma HLS pipeline II=1
				size_t l = j + 4;
				double tpre =  real[l] * cosTable[k] + imag[l] * sinTable[k];
				double tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k];
				real[l] = real[j] - tpre;
				imag[l] = imag[j] - tpim;
				real[j] += tpre;
				imag[j] += tpim;
				k += 256;
			}
		}

		// size = 16 (128 parallel)
		for (size_t i = 0; i < param_2N; i += 16) {
			#pragma HLS unroll factor=2
			int k = 0;
			for (size_t j = i; j < i + 8; j++) {
				#pragma HLS pipeline II=1
				size_t l = j + 8;
				double tpre =  real[l] * cosTable[k] + imag[l] * sinTable[k];
				double tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k];
				real[l] = real[j] - tpre;
				imag[l] = imag[j] - tpim;
				real[j] += tpre;
				imag[j] += tpim;
				k += 128;
			}
		}

		// size = 32 (64 parallel)
		for (size_t i = 0; i < param_2N; i += 32) {
			int k = 0;
			for (size_t j = i; j < i + 16; j++) {
				#pragma HLS pipeline II=1
				size_t l = j + 16;
				double tpre =  real[l] * cosTable[k] + imag[l] * sinTable[k];
				double tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k];
				real[l] = real[j] - tpre;
				imag[l] = imag[j] - tpim;
				real[j] += tpre;
				imag[j] += tpim;
				k += 64;
			}
		}

		for (size_t size = 64; size <= param_2N; size *= 2) {
			const size_t halfsize = size / 2;
			const size_t tablestep = param_2N / size;
			for (size_t i = 0; i < param_2N; i += size) {
				int k = 0;
				for (size_t j = i; j < i + halfsize; j++) {
					#pragma HLS pipeline II=1
					size_t l = j + halfsize;
					double tpre =  real[l] * cosTable[k] + imag[l] * sinTable[k];
					double tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k];
					real[l] = real[j] - tpre;
					imag[l] = imag[j] - tpim;
					real[j] += tpre;
					imag[j] += tpim;
					k += tablestep;
				}
			}
		}
	}
}