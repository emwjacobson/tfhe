#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string.h>
#include <math.h>
#include <algorithm>
#include "fpga_constants.h"

extern "C" {

	// This is a HLS implementation that models the x86-64 AVX implementation.
	void fft_transform_reverse(double *real, double *imag) {
		// Bit-reversed addressing permutation
		fft_transform_reverse_loop_1: for (uint64_t i = 0; i < param_2N; i++) {
			uint64_t j = bit_reversed[i];
			if (i < j) {
				std::swap(real[i], real[j]);
				std::swap(imag[i], imag[j]);
			}
		}

		// Size 2 merge (special)
		fft_transform_reverse_loop_2: for (uint64_t i = 0; i < param_2N; i += 2) {
			double tpre1 = real[i];
			double tpre2 = real[i + 1];
			real[i] = tpre1 + tpre2;
			real[i + 1] = tpre1 - tpre2;

			double tpim1 = imag[i];
			double tpim2 = imag[i + 1];
			imag[i] = tpim1 + tpim2;
			imag[i + 1] = tpim1 - tpim2;
		}

		// Size 4 merge (special)
		fft_transform_reverse_loop_3: for (uint64_t i = 0; i < param_2N; i += 4) {
			// Even indices
			double tpre1_even = real[i];
			double tpim1_even = imag[i];
			double tpre2_even = real[i + 2];
			double tpim2_even = imag[i + 2];
			real[i] = tpre1_even + tpre2_even;
			imag[i] = tpim1_even + tpim2_even;
			real[i + 2] = tpre1_even - tpre2_even;
			imag[i + 2] = tpim1_even - tpim2_even;


			// Odd indices
			double tpre1_odd = real[i + 1];
			double tpim1_odd = imag[i + 1];
			double tpre2_off = real[i + 3];
			double tpim2_off = imag[i + 3];
			real[i + 1] = tpre1_odd - tpim2_off;
			imag[i + 1] = tpim1_odd + tpre2_off;
			real[i + 3] = tpre1_odd + tpim2_off;
			imag[i + 3] = tpim1_odd - tpre2_off;
		}

		// Size 8 and larger merges (general)
		double *trigtables = &trig_table_reverse[0];
		fft_transform_reverse_loop_4: for (uint64_t size = 8; size <= param_2N; size *= 2) {
			const uint64_t halfsize = size >> 1;
			fft_transform_reverse_loop_5: for (uint64_t i = 0; i < param_2N; i += size) {
				uint64_t off = 0;
				fft_transform_reverse_loop_6: for (uint64_t j = 0; j < halfsize; j += 4) {
					fft_transform_reverse_loop_7: for (uint64_t k = 0; k < 4; k++) {  // To simulate x86 AVX 4-vectors
						uint64_t vi = i + j + k;  // Vector index
						uint64_t ti = off + k;    // Table index
						double re1 = real[vi + halfsize];
						double im1 = imag[vi + halfsize];
						double re2 = real[vi];
						double im2 = imag[vi];
						double tt1 = trigtables[ti];
						double tt2 = trigtables[ti + 4];
						double tpre = re1 * tt1 + im1 * tt2;
						double tpim = im1 * tt1 - re1 * tt2;
						real[vi + halfsize] = re2 - tpre;
						imag[vi + halfsize] = im2 - tpim;
						real[vi] = re2 + tpre;
						imag[vi] = im2 + tpim;
					}
					off += 8;
				}
			}
			trigtables += size;
		}

		return;
	}

}