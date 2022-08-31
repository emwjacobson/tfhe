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
			#pragma HLS dependence variable=real inter false
			#pragma HLS dependence variable=imag inter false
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

		// // size = 2
		// for (size_t i = 0; i < param_2N; i += 2) {
		// 	#pragma HLS unroll factor=64
		// 	for (size_t j = 0; j < 1; j++) {
		// 		#pragma HLS dependence variable=real inter false
		// 		#pragma HLS dependence variable=imag inter false
		// 		#pragma HLS pipeline II=1
		// 		const size_t first = j + i;
		// 		const size_t second = j + i + 1;
		// 		const int k = j * 1024;
		// 		// Load
		// 		double tprel = real[second];
		// 		double tpiml = imag[second];
		// 		double tprej = real[first];
		// 		double tpimj = imag[first];
		// 		double tpcos = cosTable[k];
		// 		double tpsin = sinTable[k];
		// 		// Calc
		// 		double calcre =  tprel * tpcos + tpiml * tpsin;
		// 		double calcim = -tprel * tpsin + tpiml * tpcos;
		// 		// Store
		// 		real[second] = tprej - calcre;
		// 		imag[second] = tpimj - calcim;
		// 		real[first] = tprej + calcre;
		// 		imag[first] = tpimj + calcim;
		// 	}
		// }

		// // size = 4
		// for (size_t i = 0; i < param_2N; i += 4) {
		// 	#pragma HLS unroll factor=64
		// 	for (size_t j = 0; j < 2; j++) {
		// 		#pragma HLS dependence variable=real inter false
		// 		#pragma HLS dependence variable=imag inter false
		// 		#pragma HLS pipeline II=1
		// 		const size_t first = j + i;
		// 		const size_t second = j + i + 2;
		// 		const int k = j * 512;
		// 		// Load
		// 		double tprel = real[second];
		// 		double tpiml = imag[second];
		// 		double tprej = real[first];
		// 		double tpimj = imag[first];
		// 		double tpcos = cosTable[k];
		// 		double tpsin = sinTable[k];
		// 		// Calc
		// 		double calcre =  tprel * tpcos + tpiml * tpsin;
		// 		double calcim = -tprel * tpsin + tpiml * tpcos;
		// 		// Store
		// 		real[second] = tprej - calcre;
		// 		imag[second] = tpimj - calcim;
		// 		real[first] = tprej + calcre;
		// 		imag[first] = tpimj + calcim;
		// 	}
		// }

		// // size = 8
		// for (size_t i = 0; i < param_2N; i += 8) {
		// 	#pragma HLS unroll factor=64
		// 	for (size_t j = 0; j < 4; j++) {
		// 		#pragma HLS dependence variable=real inter false
		// 		#pragma HLS dependence variable=imag inter false
		// 		#pragma HLS pipeline II=1
		// 		const size_t first = j + i;
		// 		const size_t second = j + i + 4;
		// 		const int k = j * 256;
		// 		// Load
		// 		double tprel = real[second];
		// 		double tpiml = imag[second];
		// 		double tprej = real[first];
		// 		double tpimj = imag[first];
		// 		double tpcos = cosTable[k];
		// 		double tpsin = sinTable[k];
		// 		// Calc
		// 		double calcre =  tprel * tpcos + tpiml * tpsin;
		// 		double calcim = -tprel * tpsin + tpiml * tpcos;
		// 		// Store
		// 		real[second] = tprej - calcre;
		// 		imag[second] = tpimj - calcim;
		// 		real[first] = tprej + calcre;
		// 		imag[first] = tpimj + calcim;
		// 	}
		// }

		// // size = 16
		// for (size_t i = 0; i < param_2N; i += 16) {
		// 	#pragma HLS unroll factor=64
		// 	for (size_t j = 0; j < 8; j++) {
		// 		#pragma HLS dependence variable=real inter false
		// 		#pragma HLS dependence variable=imag inter false
		// 		#pragma HLS pipeline II=1
		// 		const size_t first = j + i;
		// 		const size_t second = j + i + 8;
		// 		const int k = j * 128;
		// 		// Load
		// 		double tprel = real[second];
		// 		double tpiml = imag[second];
		// 		double tprej = real[first];
		// 		double tpimj = imag[first];
		// 		double tpcos = cosTable[k];
		// 		double tpsin = sinTable[k];
		// 		// Calc
		// 		double calcre =  tprel * tpcos + tpiml * tpsin;
		// 		double calcim = -tprel * tpsin + tpiml * tpcos;
		// 		// Store
		// 		real[second] = tprej - calcre;
		// 		imag[second] = tpimj - calcim;
		// 		real[first] = tprej + calcre;
		// 		imag[first] = tpimj + calcim;
		// 	}
		// }

		// // size = 32
		// for (size_t i = 0; i < param_2N; i += 32) {
		// 	#pragma HLS unroll factor=32
		// 	for (size_t j = 0; j < 16; j++) {
		// 		#pragma HLS dependence variable=real inter false
		// 		#pragma HLS dependence variable=imag inter false
		// 		#pragma HLS pipeline II=1
		// 		const size_t first = j + i;
		// 		const size_t second = j + i + 16;
		// 		const int k = j * 64;
		// 		// Load
		// 		double tprel = real[second];
		// 		double tpiml = imag[second];
		// 		double tprej = real[first];
		// 		double tpimj = imag[first];
		// 		double tpcos = cosTable[k];
		// 		double tpsin = sinTable[k];
		// 		// Calc
		// 		double calcre =  tprel * tpcos + tpiml * tpsin;
		// 		double calcim = -tprel * tpsin + tpiml * tpcos;
		// 		// Store
		// 		real[second] = tprej - calcre;
		// 		imag[second] = tpimj - calcim;
		// 		real[first] = tprej + calcre;
		// 		imag[first] = tpimj + calcim;
		// 	}
		// }

		// // size = 64
		// for (size_t i = 0; i < param_2N; i += 64) {
		// 	#pragma HLS unroll factor=16
		// 	for (size_t j = 0; j < 32; j++) {
		// 		#pragma HLS dependence variable=real inter false
		// 		#pragma HLS dependence variable=imag inter false
		// 		#pragma HLS pipeline II=1
		// 		const size_t first = j + i;
		// 		const size_t second = j + i + 32;
		// 		const int k = j * 32;
		// 		// Load
		// 		double tprel = real[second];
		// 		double tpiml = imag[second];
		// 		double tprej = real[first];
		// 		double tpimj = imag[first];
		// 		double tpcos = cosTable[k];
		// 		double tpsin = sinTable[k];
		// 		// Calc
		// 		double calcre =  tprel * tpcos + tpiml * tpsin;
		// 		double calcim = -tprel * tpsin + tpiml * tpcos;
		// 		// Store
		// 		real[second] = tprej - calcre;
		// 		imag[second] = tpimj - calcim;
		// 		real[first] = tprej + calcre;
		// 		imag[first] = tpimj + calcim;
		// 	}
		// }


		fft_size: for (size_t size = 2; size <= param_2N; size *= 2) {
			const size_t halfsize = size / 2;
			const size_t tablestep = param_2N / size;
			fft_window: for (size_t i = 0; i < param_2N; i += size) {
				fft_calc: for (size_t j = 0; j < halfsize; j++) {
					#pragma HLS dependence variable=real inter false // Amin recommended intra, inter seems to produce fewer errors?
					#pragma HLS dependence variable=imag inter false
					// #pragma HLS dependence variable=real intra false
					// #pragma HLS dependence variable=imag intra false
					#pragma HLS pipeline II=1
					const size_t first = j + i;
					const size_t second = j + i + halfsize;
					const int k = j * tablestep;
					// Load
					double tprel = real[second];
					double tpiml = imag[second];
					double tprej = real[first];
					double tpimj = imag[first];
					double tpcos = cosTable[k];
					double tpsin = sinTable[k];
					// Calc
					double calcre =  tprel * tpcos + tpiml * tpsin;
					double calcim = -tprel * tpsin + tpiml * tpcos;
					// Store
					real[second] = tprej - calcre;
					imag[second] = tpimj - calcim;
					real[first] = tprej + calcre;
					imag[first] = tpimj + calcim;
				}
			}
		}
	}
}