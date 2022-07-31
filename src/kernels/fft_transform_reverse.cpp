#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string.h>
#include <math.h>
#include "fpga_constants.h"

extern "C" {

	// This is a HLS implementation that models the x86-64 AVX implementation.
	void fft_transform_reverse(double *real, double *imag) {
		// Bit-reversed addressing permutation
		uint64_t i;
		for (i = 0; i < _2N; i++) {
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
		if (_2N >= 2) {
			for (i = 0; i < _2N; i += 2) {
				double tpre = real[i];
				double tpim = imag[i];
				real[i] += real[i + 1];
				imag[i] += imag[i + 1];
				real[i + 1] = tpre - real[i + 1];
				imag[i + 1] = tpim - imag[i + 1];
			}
		}

		// Size 4 merge (special)
		if (_2N >= 4) {
			for (i = 0; i < _2N; i += 4) {
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
		double *trigtables = trig_table_reverse;
		uint64_t size;
		for (size = 8; size <= _2N; size <<= 1) {
			uint64_t halfsize = size >> 1;
			uint64_t i;
			for (i = 0; i < _2N; i += size) {
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
			if (size == _2N)
				break;
			trigtables += size;
		}

		return;
	}

}