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
		fft_transform_reverse_loop_1: for (uint64_t i = 0; i < param_2N; i++) {
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
		if (param_2N >= 2) {
			fft_transform_reverse_loop_2: for (uint64_t i = 0; i < param_2N; i += 2) {
				double tpre = real[i];
				double tpim = imag[i];
				real[i] += real[i + 1];
				imag[i] += imag[i + 1];
				real[i + 1] = tpre - real[i + 1];
				imag[i + 1] = tpim - imag[i + 1];
			}
		}

		// Size 4 merge (special)
		if (param_2N >= 4) {
			fft_transform_reverse_loop_3: for (uint64_t i = 0; i < param_2N; i += 4) {
				// Even indices
				double tpre_even = real[i];
				double tpim_even = imag[i];
				real[i] += real[i + 2];
				imag[i] += imag[i + 2];
				real[i + 2] = tpre_even - real[i + 2];
				imag[i + 2] = tpim_even - imag[i + 2];
				// Odd indices
				double tpre_odd = real[i + 1];
				double tpim_odd = imag[i + 1];
				real[i + 1] -= imag[i + 3];
				imag[i + 1] += real[i + 3];
				tpre_odd += imag[i + 3];
				tpim_odd -= real[i + 3];
				real[i + 3] = tpre_odd;
				imag[i + 3] = tpim_odd;
			}
		}

		// Size 8 and larger merges (general)
		double *trigtables = &trig_table_reverse[0];
		fft_transform_reverse_loop_4: for (uint64_t size = 8; size <= param_2N; size *= 2) {
			// const uint64_t halfsize = size >> 1;
			fft_transform_reverse_loop_5: for (uint64_t i = 0; i < param_2N; i += size) {
				uint64_t off = 0;
				fft_transform_reverse_loop_6: for (uint64_t j = 0; j < (size >> 1); j += 4) {
					fft_transform_reverse_loop_7: for (uint64_t k = 0; k < 4; k++) {  // To simulate x86 AVX 4-vectors
						uint64_t vi = i + j + k;  // Vector index
						uint64_t ti = off + k;    // Table index
						double re = real[vi + (size >> 1)];
						double im = imag[vi + (size >> 1)];
						double tpre = re * trigtables[ti] + im * trigtables[ti + 4];
						double tpim = im * trigtables[ti] - re * trigtables[ti + 4];
						real[vi + (size >> 1)] = real[vi] - tpre;
						imag[vi + (size >> 1)] = imag[vi] - tpim;
						real[vi] += tpre;
						imag[vi] += tpim;
					}
					off += 8;
				}
			}
			trigtables += size;
		}

		return;
	}

}