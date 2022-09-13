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

	fft_mult_window: for (size_t i = 0; i < param_Ns2; i += SIZE) {
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

void load_bit_reverse(double *in_1, double *in_2, double *in_3, double *in_4, const double *in) {
	for(int i=0; i<param_2N; i++) {
		#pragma HLS pipeline II=1
		int place = i % 4;
		uint64_t j = bit_reversed[i];
		if (place == 0) {
			in_1[i] = in[j];
		} else if (place == 1) {
			in_2[i] = in[j];
		} else if (place == 2) {
			in_3[i] = in[j];
		} else if (place == 3) {
			in_4[i] = in[j];
		}
	}
}

void fft_512(double *re_out, double *im_out, const double *re_in, const double *im_in) {
	double re_temp[param_Ns2]; double im_temp[param_Ns2];

	fft_mult<2>(re_temp, im_temp, re_in, im_in);
	fft_mult<4>(re_out, im_out, re_temp, im_temp);
	fft_mult<8>(re_temp, im_temp, re_out, im_out);
	fft_mult<16>(re_out, im_out, re_temp, im_temp);
	fft_mult<32>(re_temp, im_temp, re_out, im_out);
	fft_mult<64>(re_out, im_out, re_temp, im_temp);
	fft_mult<128>(re_temp, im_temp, re_out, im_out);
	fft_mult<256>(re_out, im_out, re_temp, im_temp);
	fft_mult<512>(re_temp, im_temp, re_out, im_out);

	for(int i=0; i<param_Ns2; i++) {
		re_out[i] = re_temp[i];
		im_out[i] = im_temp[i];
	}
}

void fft_1024(double *re1_out, double *im1_out, double *re2_out, double *im2_out, const double *re1_in, const double *im1_in, const double *re2_in, const double *im2_in) {
	constexpr size_t halfsize = 1024 / 2;
	constexpr size_t tablestep = param_2N / 1024;

	fft_mult_mult: for (size_t j = 0; j < halfsize; j++) {
		#pragma HLS pipeline II=1
		const size_t first = j;
		const int k = j * tablestep;
		// Load
		double tprel = re2_in[first];
		double tpiml = im2_in[first];
		double tprej = re1_in[first];
		double tpimj = im1_in[first];
		double tpcos = cosTable[k];
		double tpsin = sinTable[k];
		// Calc
		double calcre =  tprel * tpcos + tpiml * tpsin;
		double calcim = -tprel * tpsin + tpiml * tpcos;
		// Store
		re2_out[first] = tprej - calcre;
		im2_out[first] = tpimj - calcim;
		re1_out[first] = tprej + calcre;
		im1_out[first] = tpimj + calcim;
	}
}

void fft_2048(double *re_out, double *im_out, const double *re1_in, const double *im1_in, const double *re2_in, const double *im2_in, const double *re3_in, const double *im3_in, const double *re4_in, const double *im4_in) {
	constexpr size_t halfsize = 2048 / 2;
	constexpr size_t tablestep = param_2N / 2048;

	fft_mult_mult_1: for (size_t j = 0; j < halfsize; j++) {
		#pragma HLS pipeline II=1
		const size_t first = j;
		const size_t second = j + halfsize;
		const int k = j * tablestep;
		// Load
		double tprel = re3_in[first];
		double tpiml = im3_in[first];
		double tprej = re1_in[first];
		double tpimj = im1_in[first];
		double tpcos = cosTable[k];
		double tpsin = sinTable[k];
		// Calc
		double calcre =  tprel * tpcos + tpiml * tpsin;
		double calcim = -tprel * tpsin + tpiml * tpcos;
		// Store
		re_out[second] = tprej - calcre;
		im_out[second] = tpimj - calcim;
		re_out[first] = tprej + calcre;
		im_out[first] = tpimj + calcim;
	}

	fft_mult_mult_2: for (size_t j = 0; j < halfsize; j++) {
		#pragma HLS pipeline II=1
		const size_t first = j;
		const size_t second = j + halfsize;
		const int k = j * tablestep;
		// Load
		double tprel = re4_in[first];
		double tpiml = im4_in[first];
		double tprej = re2_in[first];
		double tpimj = im2_in[first];
		double tpcos = cosTable[k];
		double tpsin = sinTable[k];
		// Calc
		double calcre =  tprel * tpcos + tpiml * tpsin;
		double calcim = -tprel * tpsin + tpiml * tpcos;
		// Store
		re_out[second+param_Ns2] = tprej - calcre;
		im_out[second+param_Ns2] = tpimj - calcim;
		re_out[first+param_Ns2] = tprej + calcre;
		im_out[first+param_Ns2] = tpimj + calcim;
	}
}

extern "C" {
	void fft_transform(double *real, double *imag) {
		double re1_1[param_Ns2]; double im1_1[param_Ns2];
		double re1_2[param_Ns2]; double im1_2[param_Ns2];
		double re1_3[param_Ns2]; double im1_3[param_Ns2];
		double re1_4[param_Ns2]; double im1_4[param_Ns2];

		double re2_1[param_Ns2]; double im2_1[param_Ns2];
		double re2_2[param_Ns2]; double im2_2[param_Ns2];
		double re2_3[param_Ns2]; double im2_3[param_Ns2];
		double re2_4[param_Ns2]; double im2_4[param_Ns2];

		#pragma HLS dataflow
		// Bit reversal
		load_bit_reverse(re1_1, re1_2, re1_3, re1_4, real);
		load_bit_reverse(im1_1, im1_2, im1_3, im1_4, imag);
		// FFT
		fft_512(re2_1, im2_1, re1_1, im1_1);
		fft_512(re2_2, im2_2, re1_2, im1_2);
		fft_512(re2_3, im2_3, re1_3, im1_3);
		fft_512(re2_4, im2_4, re1_4, im1_4);
		fft_1024(re1_1, im1_1, re1_2, im1_2, re2_1, im2_1, re2_2, im2_2);
		fft_1024(re1_3, im1_3, re1_4, im1_4, re2_3, im2_3, re2_4, im2_4);
		fft_2048(real, imag, re1_1, im1_1, re1_2, im1_2, re1_3, im1_3, re1_4, im1_4);


		// double re1[param_2N];
		// double im1[param_2N];
		// double re2[param_2N];
		// double im2[param_2N];

		// // Bit reversal
		// fft_bit_reverse: for (int i = 0; i < param_2N; i++) {
		// 	uint64_t j = bit_reversed[i];
		// 	re1[i] = real[j];
		// 	im1[i] = imag[j];
		// }

		// // At this point, temp arrays hold the bit-reversed elements
		// // Now we swap between the 2 memories to perform radix2fft
		// fft_mult<2>(re2, im2, re1, im1);
		// fft_mult<4>(re1, im1, re2, im2);
		// fft_mult<8>(re2, im2, re1, im1);
		// fft_mult<16>(re1, im1, re2, im2);
		// fft_mult<32>(re2, im2, re1, im1);
		// fft_mult<64>(re1, im1, re2, im2);
		// fft_mult<128>(re2, im2, re1, im1);
		// fft_mult<256>(re1, im1, re2, im2);
		// fft_mult<512>(re2, im2, re1, im1);
		// fft_mult<1024>(re1, im1, re2, im2);
		// fft_mult<2048>(real, imag, re1, im1);
	}
}