#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string.h>

extern "C" {
	// This is a C implementation that models the x86-64 AVX implementation.
	void fft_transform_reverse(const uint64_t _n, uint64_t *_bit_reversed, double *_trig_tables, double *_real, double *_imag) {
		// #pragma HLS INTERFACE m_axi port=n
		// #pragma HLS INTERFACE m_axi port=bit_reversed
		// #pragma HLS INTERFACE m_axi port=trig_tables
		// #pragma HLS INTERFACE m_axi port=real
		// #pragma HLS INTERFACE m_axi port=imag

		// Pull data from interface into function
		uint64_t n = _n;
		uint64_t bit_reversed[n];
		double trig_tables[n * 2];
		double real[n];
		double imag[n];

		memcpy(&bit_reversed, _bit_reversed, sizeof(uint64_t) * n);
		memcpy(&trig_tables, _trig_tables, sizeof(double) * n * 2);
		memcpy(&real, _real, sizeof(double) * n);
		memcpy(&imag, _imag, sizeof(double) * n);

		// Bit-reversed addressing permutation
		uint64_t i;
		for (i = 0; i < n; i++) {
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
		if (n >= 2) {
			for (i = 0; i < n; i += 2) {
				double tpre = real[i];
				double tpim = imag[i];
				real[i] += real[i + 1];
				imag[i] += imag[i + 1];
				real[i + 1] = tpre - real[i + 1];
				imag[i + 1] = tpim - imag[i + 1];
			}
		}

		// Size 4 merge (special)
		if (n >= 4) {
			for (i = 0; i < n; i += 4) {
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
		double *trigtables = trig_tables;
		uint64_t size;
		for (size = 8; size <= n; size <<= 1) {
			uint64_t halfsize = size >> 1;
			uint64_t i;
			for (i = 0; i < n; i += size) {
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
			if (size == n)
				break;
			trigtables += size;
		}

		memcpy(_bit_reversed, &bit_reversed, n * sizeof(uint64_t));
		memcpy(_trig_tables, &trig_tables, sizeof(double) * n * 2);
		memcpy(_real, &real, sizeof(double) * n);
		memcpy(_imag, &imag, sizeof(double) * n);

		return;
	}
}