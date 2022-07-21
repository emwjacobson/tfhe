#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string.h>
#include <math.h>

#define N 2048

extern "C" {
	double accurate_sine(uint64_t i, uint64_t n) {
		if (n % 4 != 0)
			return NAN;
		else {
			int32_t neg = 0;  // Boolean
			// Reduce to full cycle
			i %= n;
			// Reduce to half cycle
			if (i >= n / 2) {
				neg = 1;
				i -= n / 2;
			}
			// Reduce to quarter cycle
			if (i >= n / 4)
				i = n / 2 - i;
			// Reduce to eighth cycle
			double val;
			if (i * 8 < n)
				val = sin(2 * M_PI * i / n);
			else
				val = cos(2 * M_PI * (n / 4 - i) / n);
			// Apply sign
			return neg ? -val : val;
		}
	}

	void init_trig_table(double table[N*2]) {
		uint64_t k = 0;
		for (uint64_t size = 8; size <= N; size *= 2) {
			for (uint64_t i = 0; i < size / 2; i += 4) {
				uint64_t j;
				for (j = 0; j < 4; j++, k++)
					table[k] = accurate_sine(i + j + size / 4, size);  // Cosine
				for (j = 0; j < 4; j++, k++)
					table[k] = -accurate_sine(i + j, size);  // Sine
			}
			if (size == N)
				break;
		}
	}

	int32_t floor_log2(uint64_t n) {
		int32_t result = 0;
		for (; n > 1; n /= 2)
			result++;
		return result;
	}

	uint64_t reverse_bits(uint64_t x, uint32_t n) {
		uint64_t result = 0;
		uint32_t i;
		for (i = 0; i < n; i++, x >>= 1)
			result = (result << 1) | (x & 1);
		return result;
	}

	void init_bit_reversed(uint64_t bit_reversed[N]) {
		int32_t levels = floor_log2(N);
		for (uint64_t i = 0; i < N; i++)
			bit_reversed[i] = reverse_bits(i, levels);
	}

	// This is a HLS implementation that models the x86-64 AVX implementation.
	// void fft_transform_reverse(uint64_t *_bit_reversed, double *_trig_tables, double *_real, double *_imag) {
	void fft_transform_reverse(double *_real, double *_imag) {
		// #pragma HLS INTERFACE m_axi port=_bit_reversed
		// #pragma HLS INTERFACE m_axi port=_trig_tables
		#pragma HLS INTERFACE m_axi port=_real bundle=inout_real
		#pragma HLS INTERFACE m_axi port=_imag bundle=inout_imag

		// Should be implemented as ROM in final design
		static double trig_tables[N * 2];
		init_trig_table(trig_tables);

		static uint64_t bit_reversed[N];
		init_bit_reversed(bit_reversed);

		// uint64_t n = N;
		// double trig_tables[N * 2];
		// uint64_t bit_reversed[N];
		double real[N];
		double imag[N];

		// Pull data from interface into function
		// memcpy(&bit_reversed, _bit_reversed, sizeof(uint64_t) * n);
		// memcpy(&trig_tables, _trig_tables, sizeof(double) * n * 2);
		memcpy(&real, _real, sizeof(double) * N);
		memcpy(&imag, _imag, sizeof(double) * N);

		// Bit-reversed addressing permutation
		uint64_t i;
		for (i = 0; i < N; i++) {
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
			for (i = 0; i < N; i += 2) {
				double tpre = real[i];
				double tpim = imag[i];
				real[i] += real[i + 1];
				imag[i] += imag[i + 1];
				real[i + 1] = tpre - real[i + 1];
				imag[i + 1] = tpim - imag[i + 1];
			}
		}

		// Size 4 merge (special)
		if (N >= 4) {
			for (i = 0; i < N; i += 4) {
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
		for (size = 8; size <= N; size <<= 1) {
			uint64_t halfsize = size >> 1;
			uint64_t i;
			for (i = 0; i < N; i += size) {
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
			if (size == N)
				break;
			trigtables += size;
		}

		memcpy(_real, &real, sizeof(double) * N);
		memcpy(_imag, &imag, sizeof(double) * N);

		return;
	}
}