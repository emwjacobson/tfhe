#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string.h>
#include <math.h>
#include "constants.h"

extern "C" {

	// This is a HLS implementation that models the x86-64 AVX implementation.
	void fft_transform_reverse(double *_real, double *_imag) {
		#pragma HLS INTERFACE m_axi port=_real bundle=inout_real
		#pragma HLS INTERFACE m_axi port=_imag bundle=inout_imag

		double real[N];
		double imag[N];
		#pragma HLS ARRAY_PARTITION variable=real type=complete
		#pragma HLS ARRAY_PARTITION variable=imag type=complete

		// Pull data from interface into function
		memcpy(&real, _real, sizeof(double) * N);
		memcpy(&imag, _imag, sizeof(double) * N);

		// Bit-reversed addressing permutation
		for (uint64_t i = 0; i < N; i++) {
			uint64_t j = bit_reversed[i];
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
		if (N >= 2) {
			for (uint64_t i = 0; i < N; i += 2) {
				double tpre = real[i];
				double tpim = imag[i];
				real[i] = real[i] + real[i + 1];
				imag[i] = imag[i] + imag[i + 1];
				real[i + 1] = tpre - real[i + 1];
				imag[i + 1] = tpim - imag[i + 1];
			}
		}

		// Size 4 merge (special)
		if (N >= 4) {
			for (uint64_t i = 0; i < N; i += 4) {
				// Even indices
				double tpre, tpim;
				tpre = real[i];
				tpim = imag[i];
				real[i] = real[i] + real[i + 2];
				imag[i] = imag[i] + imag[i + 2];
				real[i + 2] = tpre - real[i + 2];
				imag[i + 2] = tpim - imag[i + 2];
				// Odd indices
				tpre = real[i + 1];
				tpim = imag[i + 1];
				real[i + 1] = real[i + 1] - imag[i + 3];
				imag[i + 1] = imag[i + 1] + real[i + 3];
				tpre = tpre + imag[i + 3];
				tpim = tpim - real[i + 3];
				real[i + 3] = tpre;
				imag[i + 3] = tpim;
			}
		}

		// Size 8 and larger merges (general)
		uint64_t off;
		uint64_t halfsize;
		double *trigtables = trig_tables;
		for (uint64_t size = 8; size <= N; size <<= 1) {
			for (uint64_t i = 0; i < N; i += size) {
				off = 0;
				halfsize = size >> 1;
				for (uint64_t j = 0; j < halfsize; j++) {
					uint64_t vi = i + j;  // Vector index
					uint64_t ti = off;    // Table index
					double re = real[vi + halfsize];
					double im = imag[vi + halfsize];
					double tpre = re * trigtables[ti] + im * trigtables[ti + 4];
					double tpim = im * trigtables[ti] - re * trigtables[ti + 4];
					real[vi + halfsize] = real[vi] - tpre;
					imag[vi + halfsize] = imag[vi] - tpim;
					real[vi] = real[vi] + tpre;
					imag[vi] = imag[vi] + tpim;
					off += 8;
				}
			}
			if (size == N)
				break;
			trigtables += size;
		}

		memcpy(_real, &real, sizeof(double) * N);
		memcpy(_imag, &imag, sizeof(double) * N);

		return;
	}

}