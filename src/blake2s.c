/*
 * Internal functions for BAT.
 */

/* ====================================================================== */

#include <stdint.h>
#include <string.h>

#include "blake2.h"

#include "inner.h"
#define BLAKE2_AVX2        BAT_AVX2
#define BLAKE2_LE          BAT_LE
#define BLAKE2_UNALIGNED   BAT_UNALIGNED

ALIGNED_AVX2
static const uint32_t IV[] = {
	0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
	0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19
};

#if BLAKE2_AVX2

TARGET_AVX2
static void
process_block(uint32_t *h, const uint8_t *data, uint64_t t, int f)
{
	__m128i xh0, xh1, xv0, xv1, xv2, xv3;
	__m128i xm0, xm1, xm2, xm3, xn0, xn1, xn2, xn3;
	__m128i xt0, xt1, xt2, xt3, xt4, xt5, xt6, xt7, xt8, xt9;
	__m128i xror8, xror16;

	xror8 = _mm_setr_epi8(
		1, 2, 3, 0, 5, 6, 7, 4,
		9, 10, 11, 8, 13, 14, 15, 12);
	xror16 = _mm_setr_epi8(
		2, 3, 0, 1, 6, 7, 4, 5,
		10, 11, 8, 9, 14, 15, 12, 13);

	/* Initialize state. */
	xh0 = _mm_loadu_si128((const void *)(h + 0));
	xh1 = _mm_loadu_si128((const void *)(h + 4));
	xv0 = xh0;
	xv1 = xh1;
	xv2 = _mm_loadu_si128((const void *)(IV + 0));
	xv3 = _mm_loadu_si128((const void *)(IV + 4));
	xv3 = _mm_xor_si128(xv3, _mm_setr_epi32(
		(int32_t)(uint32_t)t, (int32_t)(uint32_t)(t >> 32),
		-f, 0));

	/* Load data and move it into the proper order for the first round:
	     xm0:  0  2  4  6
	     xm1:  1  3  5  7
	     xm2:  8 10 12 14
	     xm3:  9 11 13 15 */
	xm0 = _mm_loadu_si128((const void *)(data +  0));
	xm1 = _mm_loadu_si128((const void *)(data + 16));
	xm2 = _mm_loadu_si128((const void *)(data + 32));
	xm3 = _mm_loadu_si128((const void *)(data + 48));

	xn0 = _mm_shuffle_epi32(xm0, 0xD8);
	xn1 = _mm_shuffle_epi32(xm1, 0xD8);
	xm0 = _mm_unpacklo_epi64(xn0, xn1);
	xm1 = _mm_unpackhi_epi64(xn0, xn1);

	xn2 = _mm_shuffle_epi32(xm2, 0xD8);
	xn3 = _mm_shuffle_epi32(xm3, 0xD8);
	xm2 = _mm_unpacklo_epi64(xn2, xn3);
	xm3 = _mm_unpackhi_epi64(xn2, xn3);

#define G4(xx, xy)   do { \
		__m128i xtg; \
		xv0 = _mm_add_epi32(xv0, _mm_add_epi32(xv1, xx)); \
		xv3 = _mm_shuffle_epi8(_mm_xor_si128(xv0, xv3), xror16); \
		xv2 = _mm_add_epi32(xv2, xv3); \
		xtg = _mm_xor_si128(xv1, xv2); \
		xv1 = _mm_or_si128( \
			_mm_srli_epi32(xtg, 12), _mm_slli_epi32(xtg, 20)); \
		xv0 = _mm_add_epi32(xv0, _mm_add_epi32(xv1, xy)); \
		xv3 = _mm_shuffle_epi8(_mm_xor_si128(xv0, xv3), xror8); \
		xv2 = _mm_add_epi32(xv2, xv3); \
		xtg = _mm_xor_si128(xv1, xv2); \
		xv1 = _mm_or_si128( \
			_mm_srli_epi32(xtg, 7), _mm_slli_epi32(xtg, 25)); \
	} while (0)

#define ROUND(i0, i1, i2, i3)   do { \
		G4(i0, i1); \
		xv1 = _mm_shuffle_epi32(xv1, 0x39); \
		xv2 = _mm_shuffle_epi32(xv2, 0x4E); \
		xv3 = _mm_shuffle_epi32(xv3, 0x93); \
		G4(i2, i3); \
		xv1 = _mm_shuffle_epi32(xv1, 0x93); \
		xv2 = _mm_shuffle_epi32(xv2, 0x4E); \
		xv3 = _mm_shuffle_epi32(xv3, 0x39); \
	} while (0)

	/* round 0 */
	ROUND(xm0, xm1, xm2, xm3);

	/* round 1 */
	xt0 = _mm_shuffle_epi32(xm0, 0x00);
	xt1 = _mm_shuffle_epi32(xm0, 0xC8);
	xt2 = _mm_shuffle_epi32(xm1, 0x70);
	xt3 = _mm_shuffle_epi32(xm1, 0x80);
	xt4 = _mm_shuffle_epi32(xm2, 0x01);
	xt5 = _mm_shuffle_epi32(xm2, 0x02);
	xt6 = _mm_shuffle_epi32(xm2, 0x03);
	xt7 = _mm_shuffle_epi32(xm3, 0x80);
	xt8 = _mm_shuffle_epi32(xm3, 0x10);
	xt9 = _mm_shuffle_epi32(xm3, 0x30);
	xn0 = _mm_blend_epi32(
		_mm_blend_epi32(xt6, xt1, 0x02),
		xt7, 0x0C);
	xn1 = _mm_blend_epi32(
		_mm_blend_epi32(xt4, xt9, 0x04),
		xt1, 0x08);
	xn2 = _mm_blend_epi32(
		_mm_blend_epi32(xt3, xt0, 0x02),
		xt8, 0x04);
	xn3 = _mm_blend_epi32(
		_mm_blend_epi32(xt5, xm0, 0x02),
		xt2, 0x0C);
	ROUND(xn0, xn1, xn2, xn3);

	/* round 2 */
	xt0 = _mm_shuffle_epi32(xn0, 0x40);
	xt1 = _mm_shuffle_epi32(xn0, 0x80);
	xt2 = _mm_shuffle_epi32(xn1, 0x80);
	xt3 = _mm_shuffle_epi32(xn1, 0x0D);
	xt4 = _mm_shuffle_epi32(xn2, 0x04);
	xt5 = _mm_shuffle_epi32(xn2, 0x32);
	xt6 = _mm_shuffle_epi32(xn3, 0x10);
	xt7 = _mm_shuffle_epi32(xn3, 0x2C);
	xm0 = _mm_blend_epi32(
		_mm_blend_epi32(xt5, xt6, 0x02),
		xt2, 0x08);
	xm1 = _mm_blend_epi32(
		_mm_blend_epi32(xt3, xt4, 0x02),
		_mm_blend_epi32(xt6, xn0, 0x08), 0x0C);
	xm2 = _mm_blend_epi32(
		_mm_blend_epi32(xt2, xt7, 0x06),
		xt1, 0x08);
	xm3 = _mm_blend_epi32(
		_mm_blend_epi32(xt0, xt3, 0x02),
		xt4, 0x04);
	ROUND(xm0, xm1, xm2, xm3);

	/* round 3 */
	xt0 = _mm_shuffle_epi32(xm0, 0x10);
	xt1 = _mm_shuffle_epi32(xm0, 0xC8);
	xt2 = _mm_shuffle_epi32(xm1, 0x10);
	xt3 = _mm_shuffle_epi32(xm1, 0x32);
	xt4 = _mm_shuffle_epi32(xm2, 0x03);
	xt5 = _mm_shuffle_epi32(xm2, 0x06);
	xt6 = _mm_shuffle_epi32(xm3, 0x39);
	xn0 = _mm_blend_epi32(
		_mm_blend_epi32(xt5, xt3, 0x04),
		xt0, 0x08);
	xn1 = _mm_blend_epi32(
		_mm_blend_epi32(xt4, xt6, 0x0A),
		xt0, 0x04);
	xn2 = _mm_blend_epi32(
		_mm_blend_epi32(xt3, xt1, 0x0A),
		xt6, 0x04);
	xn3 = _mm_blend_epi32(
		_mm_blend_epi32(xt6, xt4, 0x02),
		xt2, 0x0C);
	ROUND(xn0, xn1, xn2, xn3);

	/* round 4 */
	xt0 = _mm_shuffle_epi32(xn0, 0x80);
	xt1 = _mm_shuffle_epi32(xn0, 0x4C);
	xt2 = _mm_shuffle_epi32(xn1, 0x09);
	xt3 = _mm_shuffle_epi32(xn1, 0x03);
	xt4 = _mm_shuffle_epi32(xn2, 0x04);
	xt5 = _mm_shuffle_epi32(xn3, 0x40);
	xt6 = _mm_shuffle_epi32(xn3, 0x32);
	xm0 = _mm_blend_epi32(
		_mm_blend_epi32(xn1, xt4, 0x06),
		xt5, 0x08);
	xm1 = _mm_blend_epi32(
		_mm_blend_epi32(xt6, xt0, 0x02),
		xn2, 0x0C);
	xm2 = _mm_blend_epi32(
		_mm_blend_epi32(xt3, xt1, 0x0A),
		xt5, 0x04);
	xm3 = _mm_blend_epi32(
		_mm_blend_epi32(xt2, xt6, 0x04),
		xt0, 0x08);
	ROUND(xm0, xm1, xm2, xm3);

	/* round 5 */
	xt0 = _mm_shuffle_epi32(xm0, 0x04);
	xt1 = _mm_shuffle_epi32(xm0, 0x0E);
	xt2 = _mm_shuffle_epi32(xm1, 0x04);
	xt3 = _mm_shuffle_epi32(xm1, 0x32);
	xt4 = _mm_shuffle_epi32(xm2, 0x08);
	xt5 = _mm_shuffle_epi32(xm2, 0xD0);
	xt6 = _mm_shuffle_epi32(xm3, 0x01);
	xt7 = _mm_shuffle_epi32(xm3, 0x83);
	xn0 = _mm_blend_epi32(
		_mm_blend_epi32(xt1, xt4, 0x02),
		_mm_blend_epi32(xt2, xt7, 0x08), 0x0C);
	xn1 = _mm_blend_epi32(
		_mm_blend_epi32(xt6, xt1, 0x02),
		xt5, 0x0C);
	xn2 = _mm_blend_epi32(
		_mm_blend_epi32(xt3, xt2, 0x02),
		xt6, 0x08);
	xn3 = _mm_blend_epi32(
		_mm_blend_epi32(xt7, xt0, 0x0A),
		xt4, 0x04);
	ROUND(xn0, xn1, xn2, xn3);

	/* round 6 */
	xt0 = _mm_shuffle_epi32(xn0, 0xC6);
	xt1 = _mm_shuffle_epi32(xn1, 0x40);
	xt2 = _mm_shuffle_epi32(xn1, 0x8C);
	xt3 = _mm_shuffle_epi32(xn2, 0x09);
	xt4 = _mm_shuffle_epi32(xn2, 0x0C);
	xt5 = _mm_shuffle_epi32(xn3, 0x01);
	xt6 = _mm_shuffle_epi32(xn3, 0x30);
	xm0 = _mm_blend_epi32(
		_mm_blend_epi32(xt1, xt4, 0x0A),
		xn3, 0x04);
	xm1 = _mm_blend_epi32(
		_mm_blend_epi32(xt5, xt3, 0x02),
		xt1, 0x08);
	xm2 = _mm_blend_epi32(xt0, xt6, 0x04);
	xm3 = _mm_blend_epi32(
		_mm_blend_epi32(xt3, xt2, 0x0A),
		xt0, 0x04);
	ROUND(xm0, xm1, xm2, xm3);

	/* round 7 */
	xt0 = _mm_shuffle_epi32(xm0, 0x0C);
	xt1 = _mm_shuffle_epi32(xm0, 0x18);
	xt2 = _mm_shuffle_epi32(xm1, 0xC2);
	xt3 = _mm_shuffle_epi32(xm2, 0x10);
	xt4 = _mm_shuffle_epi32(xm2, 0xB0);
	xt5 = _mm_shuffle_epi32(xm3, 0x40);
	xt6 = _mm_shuffle_epi32(xm3, 0x83);
	xn0 = _mm_blend_epi32(
		_mm_blend_epi32(xt2, xt5, 0x0A),
		xt0, 0x04);
	xn1 = _mm_blend_epi32(
		_mm_blend_epi32(xt6, xt1, 0x06),
		xt4, 0x08);
	xn2 = _mm_blend_epi32(
		_mm_blend_epi32(xm1, xt4, 0x04),
		xt6, 0x08);
	xn3 = _mm_blend_epi32(
		_mm_blend_epi32(xt3, xt0, 0x02),
		xt2, 0x08);
	ROUND(xn0, xn1, xn2, xn3);

	/* round 8 */
	xt0 = _mm_shuffle_epi32(xn0, 0x02);
	xt1 = _mm_shuffle_epi32(xn0, 0x34);
	xt2 = _mm_shuffle_epi32(xn1, 0x0C);
	xt3 = _mm_shuffle_epi32(xn2, 0x03);
	xt4 = _mm_shuffle_epi32(xn2, 0x81);
	xt5 = _mm_shuffle_epi32(xn3, 0x02);
	xt6 = _mm_shuffle_epi32(xn3, 0xD0);
	xm0 = _mm_blend_epi32(
		_mm_blend_epi32(xt5, xn1, 0x02),
		xt2, 0x04);
	xm1 = _mm_blend_epi32(
		_mm_blend_epi32(xt4, xt2, 0x02),
		xt1, 0x04);
	xm2 = _mm_blend_epi32(
		_mm_blend_epi32(xt0, xn1, 0x04),
		xt6, 0x08);
	xm3 = _mm_blend_epi32(
		_mm_blend_epi32(xt3, xt1, 0x02),
		xt6, 0x04);
	ROUND(xm0, xm1, xm2, xm3);

	/* round 9 */
	xt0 = _mm_shuffle_epi32(xm0, 0xC6);
	xt1 = _mm_shuffle_epi32(xm1, 0x2C);
	xt2 = _mm_shuffle_epi32(xm2, 0x40);
	xt3 = _mm_shuffle_epi32(xm2, 0x83);
	xt4 = _mm_shuffle_epi32(xm3, 0xD8);
	xn0 = _mm_blend_epi32(
		_mm_blend_epi32(xt3, xt1, 0x02),
		xt4, 0x04);
	xn1 = _mm_blend_epi32(xt4, xt0, 0x04);
	xn2 = _mm_blend_epi32(
		_mm_blend_epi32(xm1, xt1, 0x04),
		xt2, 0x08);
	xn3 = _mm_blend_epi32(xt0, xt2, 0x04);
	ROUND(xn0, xn1, xn2, xn3);

#undef G4
#undef ROUND

	xh0 = _mm_xor_si128(xh0, _mm_xor_si128(xv0, xv2));
	xh1 = _mm_xor_si128(xh1, _mm_xor_si128(xv1, xv3));
	_mm_storeu_si128((void *)(h + 0), xh0);
	_mm_storeu_si128((void *)(h + 4), xh1);
}

/*
 * Optimized AVX2 implementation for blake2s_expand().
 *
 * Input buffer data_x2[] has size 128 bytes; it is filled with two
 * interlaced instances of the label, initial block counter and seed in
 * their proper positions, and zero elsewhere.
 *
 * In data_x2[], the counter fields are supposed to be already set for the
 * two first blocks.
 *
 * 'data_len' is the message length (16 + length of seed, not the length
 * of the duplicated buffer data_x2).
 *
 * This function produces dst_len bytes. dst_len MUST be non-zero,
 * and a multiple of 64 bytes.
 *
 * The function internally increments the block counters over 32 bits only,
 * without carry propagation. The caller is responsible for calling this
 * function with initial counter and dst_len values that ensure that no
 * carry propagation is missed.
 */
TARGET_AVX2
static void
expand_inner_x2(uint8_t *dst, size_t dst_len,
	const uint8_t *data_x2, size_t data_len)
{
	/* Initial value, duplicated for AVX2 parallelism. */
	ALIGNED_AVX2
	static const uint32_t IV_x2[] = {
		0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
		0x6A09E667, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
		0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19,
		0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19
	};

	/* Initial state, duplicated for AVX2 parallelism, with a
	   personalization block for output 32 bytes. */
	ALIGNED_AVX2
	static const uint32_t hinit_out32_x2[] = {
		0x6A09E667 ^ 0x01010020, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
		0x6A09E667 ^ 0x01010020, 0xBB67AE85, 0x3C6EF372, 0xA54FF53A,
		0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19,
		0x510E527F, 0x9B05688C, 0x1F83D9AB, 0x5BE0CD19
	};

	__m256i xh0, xh1, xv0, xv1, xv2, xv3;
	__m256i xm0, xm1, xm2, xm3, xn0, xn1, xn2, xn3;
	__m256i xt0, xt1, xt2, xt3, xt4, xt5, xt6, xt7, xt8, xt9;
	__m256i xror8, xror16, xca2;

	xror8 = _mm256_setr_epi8(
		1, 2, 3, 0, 5, 6, 7, 4,
		9, 10, 11, 8, 13, 14, 15, 12,
		1, 2, 3, 0, 5, 6, 7, 4,
		9, 10, 11, 8, 13, 14, 15, 12);
	xror16 = _mm256_setr_epi8(
		2, 3, 0, 1, 6, 7, 4, 5,
		10, 11, 8, 9, 14, 15, 12, 13,
		2, 3, 0, 1, 6, 7, 4, 5,
		10, 11, 8, 9, 14, 15, 12, 13);
	xca2 = _mm256_setr_epi32(0, 2, 0, 0, 0, 2, 0, 0);

	/* Initialize state. */
	xh0 = _mm256_loadu_si256((const void *)(hinit_out32_x2 + 0));
	xh1 = _mm256_loadu_si256((const void *)(hinit_out32_x2 + 8));

	/* Load data and move it into the proper order for the first round:
	     xm0:  0  2  4  6
	     xm1:  1  3  5  7
	     xm2:  8 10 12 14
	     xm3:  9 11 13 15 */
	xm0 = _mm256_loadu_si256((const void *)(data_x2 +  0));
	xm1 = _mm256_loadu_si256((const void *)(data_x2 + 32));
	xm2 = _mm256_loadu_si256((const void *)(data_x2 + 64));
	xm3 = _mm256_loadu_si256((const void *)(data_x2 + 96));

	xn0 = _mm256_shuffle_epi32(xm0, 0xD8);
	xn1 = _mm256_shuffle_epi32(xm1, 0xD8);
	xm0 = _mm256_unpacklo_epi64(xn0, xn1);
	xm1 = _mm256_unpackhi_epi64(xn0, xn1);

	xn2 = _mm256_shuffle_epi32(xm2, 0xD8);
	xn3 = _mm256_shuffle_epi32(xm3, 0xD8);
	xm2 = _mm256_unpacklo_epi64(xn2, xn3);
	xm3 = _mm256_unpackhi_epi64(xn2, xn3);

	for (;;) {
		/* Each loop iteration computes two BLAKE2s in parallel,
		   in the low and high lanes, respectively. */

		/* Initialize round state. */
		xv0 = xh0;
		xv1 = xh1;
		xv2 = _mm256_loadu_si256((const void *)(IV_x2 + 0));
		xv3 = _mm256_loadu_si256((const void *)(IV_x2 + 8));
		xv3 = _mm256_xor_si256(xv3,
			_mm256_setr_epi64x(
				data_len, 0xFFFFFFFF,
				data_len, 0xFFFFFFFF));

#define G4(xx, xy)   do { \
		__m256i xtg; \
		xv0 = _mm256_add_epi32(xv0, _mm256_add_epi32(xv1, xx)); \
		xv3 = _mm256_shuffle_epi8(_mm256_xor_si256(xv0, xv3), xror16); \
		xv2 = _mm256_add_epi32(xv2, xv3); \
		xtg = _mm256_xor_si256(xv1, xv2); \
		xv1 = _mm256_or_si256( \
			_mm256_srli_epi32(xtg, 12), \
			_mm256_slli_epi32(xtg, 20)); \
		xv0 = _mm256_add_epi32(xv0, _mm256_add_epi32(xv1, xy)); \
		xv3 = _mm256_shuffle_epi8(_mm256_xor_si256(xv0, xv3), xror8); \
		xv2 = _mm256_add_epi32(xv2, xv3); \
		xtg = _mm256_xor_si256(xv1, xv2); \
		xv1 = _mm256_or_si256( \
			_mm256_srli_epi32(xtg, 7), \
			_mm256_slli_epi32(xtg, 25)); \
	} while (0)

#define ROUND(i0, i1, i2, i3)   do { \
		G4(i0, i1); \
		xv1 = _mm256_shuffle_epi32(xv1, 0x39); \
		xv2 = _mm256_shuffle_epi32(xv2, 0x4E); \
		xv3 = _mm256_shuffle_epi32(xv3, 0x93); \
		G4(i2, i3); \
		xv1 = _mm256_shuffle_epi32(xv1, 0x93); \
		xv2 = _mm256_shuffle_epi32(xv2, 0x4E); \
		xv3 = _mm256_shuffle_epi32(xv3, 0x39); \
	} while (0)

		/* round 0 */
		ROUND(xm0, xm1, xm2, xm3);

		/* round 1 */
		xt0 = _mm256_shuffle_epi32(xm0, 0x00);
		xt1 = _mm256_shuffle_epi32(xm0, 0xC8);
		xt2 = _mm256_shuffle_epi32(xm1, 0x70);
		xt3 = _mm256_shuffle_epi32(xm1, 0x80);
		xt4 = _mm256_shuffle_epi32(xm2, 0x01);
		xt5 = _mm256_shuffle_epi32(xm2, 0x02);
		xt6 = _mm256_shuffle_epi32(xm2, 0x03);
		xt7 = _mm256_shuffle_epi32(xm3, 0x80);
		xt8 = _mm256_shuffle_epi32(xm3, 0x10);
		xt9 = _mm256_shuffle_epi32(xm3, 0x30);
		xn0 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt6, xt1, 0x22),
			xt7, 0xCC);
		xn1 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt4, xt9, 0x44),
			xt1, 0x88);
		xn2 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt3, xt0, 0x22),
			xt8, 0x44);
		xn3 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt5, xm0, 0x22),
			xt2, 0xCC);
		ROUND(xn0, xn1, xn2, xn3);

		/* round 2 */
		xt0 = _mm256_shuffle_epi32(xn0, 0x40);
		xt1 = _mm256_shuffle_epi32(xn0, 0x80);
		xt2 = _mm256_shuffle_epi32(xn1, 0x80);
		xt3 = _mm256_shuffle_epi32(xn1, 0x0D);
		xt4 = _mm256_shuffle_epi32(xn2, 0x04);
		xt5 = _mm256_shuffle_epi32(xn2, 0x32);
		xt6 = _mm256_shuffle_epi32(xn3, 0x10);
		xt7 = _mm256_shuffle_epi32(xn3, 0x2C);
		xm0 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt5, xt6, 0x22),
			xt2, 0x88);
		xm1 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt3, xt4, 0x22),
			_mm256_blend_epi32(xt6, xn0, 0x88), 0xCC);
		xm2 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt2, xt7, 0x66),
			xt1, 0x88);
		xm3 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt0, xt3, 0x22),
			xt4, 0x44);
		ROUND(xm0, xm1, xm2, xm3);

		/* round 3 */
		xt0 = _mm256_shuffle_epi32(xm0, 0x10);
		xt1 = _mm256_shuffle_epi32(xm0, 0xC8);
		xt2 = _mm256_shuffle_epi32(xm1, 0x10);
		xt3 = _mm256_shuffle_epi32(xm1, 0x32);
		xt4 = _mm256_shuffle_epi32(xm2, 0x03);
		xt5 = _mm256_shuffle_epi32(xm2, 0x06);
		xt6 = _mm256_shuffle_epi32(xm3, 0x39);
		xn0 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt5, xt3, 0x44),
			xt0, 0x88);
		xn1 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt4, xt6, 0xAA),
			xt0, 0x44);
		xn2 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt3, xt1, 0xAA),
			xt6, 0x44);
		xn3 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt6, xt4, 0x22),
			xt2, 0xCC);
		ROUND(xn0, xn1, xn2, xn3);

		/* round 4 */
		xt0 = _mm256_shuffle_epi32(xn0, 0x80);
		xt1 = _mm256_shuffle_epi32(xn0, 0x4C);
		xt2 = _mm256_shuffle_epi32(xn1, 0x09);
		xt3 = _mm256_shuffle_epi32(xn1, 0x03);
		xt4 = _mm256_shuffle_epi32(xn2, 0x04);
		xt5 = _mm256_shuffle_epi32(xn3, 0x40);
		xt6 = _mm256_shuffle_epi32(xn3, 0x32);
		xm0 = _mm256_blend_epi32(
			_mm256_blend_epi32(xn1, xt4, 0x66),
			xt5, 0x88);
		xm1 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt6, xt0, 0x22),
			xn2, 0xCC);
		xm2 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt3, xt1, 0xAA),
			xt5, 0x44);
		xm3 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt2, xt6, 0x44),
			xt0, 0x88);
		ROUND(xm0, xm1, xm2, xm3);

		/* round 5 */
		xt0 = _mm256_shuffle_epi32(xm0, 0x04);
		xt1 = _mm256_shuffle_epi32(xm0, 0x0E);
		xt2 = _mm256_shuffle_epi32(xm1, 0x04);
		xt3 = _mm256_shuffle_epi32(xm1, 0x32);
		xt4 = _mm256_shuffle_epi32(xm2, 0x08);
		xt5 = _mm256_shuffle_epi32(xm2, 0xD0);
		xt6 = _mm256_shuffle_epi32(xm3, 0x01);
		xt7 = _mm256_shuffle_epi32(xm3, 0x83);
		xn0 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt1, xt4, 0x22),
			_mm256_blend_epi32(xt2, xt7, 0x88), 0xCC);
		xn1 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt6, xt1, 0x22),
			xt5, 0xCC);
		xn2 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt3, xt2, 0x22),
			xt6, 0x88);
		xn3 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt7, xt0, 0xAA),
			xt4, 0x44);
		ROUND(xn0, xn1, xn2, xn3);

		/* round 6 */
		xt0 = _mm256_shuffle_epi32(xn0, 0xC6);
		xt1 = _mm256_shuffle_epi32(xn1, 0x40);
		xt2 = _mm256_shuffle_epi32(xn1, 0x8C);
		xt3 = _mm256_shuffle_epi32(xn2, 0x09);
		xt4 = _mm256_shuffle_epi32(xn2, 0x0C);
		xt5 = _mm256_shuffle_epi32(xn3, 0x01);
		xt6 = _mm256_shuffle_epi32(xn3, 0x30);
		xm0 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt1, xt4, 0xAA),
			xn3, 0x44);
		xm1 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt5, xt3, 0x22),
			xt1, 0x88);
		xm2 = _mm256_blend_epi32(xt0, xt6, 0x44);
		xm3 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt3, xt2, 0xAA),
			xt0, 0x44);
		ROUND(xm0, xm1, xm2, xm3);

		/* round 7 */
		xt0 = _mm256_shuffle_epi32(xm0, 0x0C);
		xt1 = _mm256_shuffle_epi32(xm0, 0x18);
		xt2 = _mm256_shuffle_epi32(xm1, 0xC2);
		xt3 = _mm256_shuffle_epi32(xm2, 0x10);
		xt4 = _mm256_shuffle_epi32(xm2, 0xB0);
		xt5 = _mm256_shuffle_epi32(xm3, 0x40);
		xt6 = _mm256_shuffle_epi32(xm3, 0x83);
		xn0 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt2, xt5, 0xAA),
			xt0, 0x44);
		xn1 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt6, xt1, 0x66),
			xt4, 0x88);
		xn2 = _mm256_blend_epi32(
			_mm256_blend_epi32(xm1, xt4, 0x44),
			xt6, 0x88);
		xn3 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt3, xt0, 0x22),
			xt2, 0x88);
		ROUND(xn0, xn1, xn2, xn3);

		/* round 8 */
		xt0 = _mm256_shuffle_epi32(xn0, 0x02);
		xt1 = _mm256_shuffle_epi32(xn0, 0x34);
		xt2 = _mm256_shuffle_epi32(xn1, 0x0C);
		xt3 = _mm256_shuffle_epi32(xn2, 0x03);
		xt4 = _mm256_shuffle_epi32(xn2, 0x81);
		xt5 = _mm256_shuffle_epi32(xn3, 0x02);
		xt6 = _mm256_shuffle_epi32(xn3, 0xD0);
		xm0 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt5, xn1, 0x22),
			xt2, 0x44);
		xm1 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt4, xt2, 0x22),
			xt1, 0x44);
		xm2 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt0, xn1, 0x44),
			xt6, 0x88);
		xm3 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt3, xt1, 0x22),
			xt6, 0x44);
		ROUND(xm0, xm1, xm2, xm3);

		/* round 9 */
		xt0 = _mm256_shuffle_epi32(xm0, 0xC6);
		xt1 = _mm256_shuffle_epi32(xm1, 0x2C);
		xt2 = _mm256_shuffle_epi32(xm2, 0x40);
		xt3 = _mm256_shuffle_epi32(xm2, 0x83);
		xt4 = _mm256_shuffle_epi32(xm3, 0xD8);
		xn0 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt3, xt1, 0x22),
			xt4, 0x44);
		xn1 = _mm256_blend_epi32(xt4, xt0, 0x44);
		xn2 = _mm256_blend_epi32(
			_mm256_blend_epi32(xm1, xt1, 0x44),
			xt2, 0x88);
		xn3 = _mm256_blend_epi32(xt0, xt2, 0x44);
		ROUND(xn0, xn1, xn2, xn3);
#undef G4
#undef ROUND

		/* Finalize computation and store output. The output must
		   be deinterlaced since output blocks are supposed to
		   be consecutive . */
		xt0 = _mm256_xor_si256(xh0, _mm256_xor_si256(xv0, xv2));
		xt1 = _mm256_xor_si256(xh1, _mm256_xor_si256(xv1, xv3));
		xt2 = _mm256_permute2x128_si256(xt0, xt1, 0x20);
		xt3 = _mm256_permute2x128_si256(xt0, xt1, 0x31);
		_mm256_storeu_si256((void *)(dst +  0), xt2);
		_mm256_storeu_si256((void *)(dst + 32), xt3);

		dst += 64;
		dst_len -= 64;
		if (dst_len == 0) {
			break;
		}

		/* Put back message words in initial order */
		xt0 = _mm256_shuffle_epi32(xn0, 0x01);
		xt1 = _mm256_shuffle_epi32(xn0, 0x83);
		xt2 = _mm256_shuffle_epi32(xn1, 0x10);
		xt3 = _mm256_shuffle_epi32(xn1, 0xB0);
		xt4 = _mm256_shuffle_epi32(xn2, 0x39);
		xt5 = _mm256_shuffle_epi32(xn3, 0x63);
		xm0 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt5, xt2, 0x66),
			xt3, 0x88);
		xm1 = _mm256_blend_epi32(
			_mm256_blend_epi32(xt1, xt4, 0x22),
			xt3, 0x44);
		xm2 = _mm256_blend_epi32(xt0, xt5, 0xCC);
		xm3 = _mm256_blend_epi32(xt4, xt5, 0x22);

		/* Increment block counter in the message.
		   Nominally, the counter is 64 bits, but we only
		   increment the low 32 bits; the caller is responsible
		   for setting the high half and calling us with values
		   that won't overflow. */
		xm0 = _mm256_add_epi32(xm0, xca2);
	}
}

#else

static void
process_block(uint32_t *h, const uint8_t *data, uint64_t t, int f)
{
	uint32_t v[16], m[16];
	int i;

	memcpy(v, h, 8 * sizeof(uint32_t));
	memcpy(v + 8, IV, sizeof IV);
	v[12] ^= (uint32_t)t;
	v[13] ^= (uint32_t)(t >> 32);
	if (f) {
		v[14] = ~v[14];
	}

#if BLAKE2_LE
	memcpy(m, data, sizeof m);
#else
	for (i = 0; i < 16; i ++) {
		m[i] = dec32le(data + (i << 2));
	}
#endif

#define ROR(x, n)   (((x) << (32 - (n))) | ((x) >> (n)))

#define G(a, b, c, d, x, y)   do { \
		v[a] += v[b] + (x); \
		v[d] = ROR(v[d] ^ v[a], 16); \
		v[c] += v[d]; \
		v[b] = ROR(v[b] ^ v[c], 12); \
		v[a] += v[b] + (y); \
		v[d] = ROR(v[d] ^ v[a], 8); \
		v[c] += v[d]; \
		v[b] = ROR(v[b] ^ v[c], 7); \
	} while (0)

#define ROUND(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, sA, sB, sC, sD, sE, sF) \
	do { \
		G(0, 4,  8, 12, m[s0], m[s1]); \
		G(1, 5,  9, 13, m[s2], m[s3]); \
		G(2, 6, 10, 14, m[s4], m[s5]); \
		G(3, 7, 11, 15, m[s6], m[s7]); \
		G(0, 5, 10, 15, m[s8], m[s9]); \
		G(1, 6, 11, 12, m[sA], m[sB]); \
		G(2, 7,  8, 13, m[sC], m[sD]); \
		G(3, 4,  9, 14, m[sE], m[sF]); \
	} while (0)

	ROUND( 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15);
	ROUND(14, 10,  4,  8,  9, 15, 13,  6,  1, 12,  0,  2, 11,  7,  5,  3);
	ROUND(11,  8, 12,  0,  5,  2, 15, 13, 10, 14,  3,  6,  7,  1,  9,  4);
	ROUND( 7,  9,  3,  1, 13, 12, 11, 14,  2,  6,  5, 10,  4,  0, 15,  8);
	ROUND( 9,  0,  5,  7,  2,  4, 10, 15, 14,  1, 11, 12,  6,  8,  3, 13);
	ROUND( 2, 12,  6, 10,  0, 11,  8,  3,  4, 13,  7,  5, 15, 14,  1,  9);
	ROUND(12,  5,  1, 15, 14, 13,  4, 10,  0,  7,  6,  3,  9,  2,  8, 11);
	ROUND(13, 11,  7, 14, 12,  1,  3,  9,  5,  0, 15,  4,  8,  6,  2, 10);
	ROUND( 6, 15, 14,  9, 11,  3,  0,  8, 12,  2, 13,  7,  1,  4, 10,  5);
	ROUND(10,  2,  8,  4,  7,  6,  1,  5, 15, 11,  9, 14,  3, 12, 13,  0);

#undef ROR
#undef G
#undef ROUND

	for (i = 0; i < 8; i ++) {
		h[i] ^= v[i] ^ v[i + 8];
	}
}

#endif

/*
 * State rules:
 *
 *   buf    buffered data
 *   h      current state
 *   ctr    number of bytes injected so far
 *
 * Initially, ctr == 0 and h contains the XOR of IV and parameter block;
 * buf[] is empty. For any ctr > 0, buf[] is non-empty; it might contain
 * a full block worth of data (processing of the block is delayed until
 * we know whether this is the final block or not).
 *
 * If a key is injected, then it counts as a first full block.
 */

/* see blake2.h */
void
blake2s_init(blake2s_context *bc, size_t out_len)
{
	memcpy(bc->h, IV, sizeof bc->h);
	bc->h[0] ^= 0x01010000 ^ (uint32_t)out_len;
	bc->ctr = 0;
	bc->out_len = out_len;
}

/* see blake2.h */
void
blake2s_init_key(blake2s_context *bc, size_t out_len,
	const void *key, size_t key_len)
{
	blake2s_init(bc, out_len);
	if (key_len > 0) {
		bc->h[0] ^= (uint32_t)key_len << 8;
		memcpy(bc->buf, key, key_len);
		memset(bc->buf + key_len, 0, (sizeof bc->buf) - key_len);
		bc->ctr = sizeof bc->buf;
	}
}

/* see blake2.h */
void
blake2s_update(blake2s_context *bc, const void *data, size_t len)
{
	uint64_t ctr;
	size_t p;

	/* Special case: if no input data, return immediately. */
	if (len == 0) {
		return;
	}

	ctr = bc->ctr;

	/* First complete the current block, if not already full. */
	p = (size_t)ctr & ((sizeof bc->buf) - 1);
	if (ctr == 0 || p != 0) {
		/* buffer is not full */
		size_t clen;

		clen = sizeof bc->buf - p;
		if (clen >= len) {
			memcpy(bc->buf + p, data, len);
			bc->ctr = ctr + len;
			return;
		}
		memcpy(bc->buf + p, data, clen);
		ctr += clen;
		data = (const uint8_t *)data + clen;
		len -= clen;
	}

	/* Process the buffered block. */
	process_block(bc->h, bc->buf, ctr, 0);

	/* Process all subsequent full blocks, except the last. */
	while (len > sizeof bc->buf) {
		ctr += sizeof bc->buf;
		process_block(bc->h, data, ctr, 0);
		data = (const uint8_t *)data + sizeof bc->buf;
		len -= sizeof bc->buf;
	}

	/* Copy the last block (possibly partial) into the buffer. */
	memcpy(bc->buf, data, len);
	bc->ctr = ctr + len;
}

/* see blake2.h */
void
blake2s_final(blake2s_context *bc, void *dst)
{
#if !BLAKE2_LE
	int i;
	uint8_t tmp[32];
#endif
	size_t p;

	/* Pad the current block with zeros, if not full. If the
	   buffer is empty (no key, no data) then fill it with zeros
	   as well. */
	p = (size_t)bc->ctr & ((sizeof bc->buf) - 1);
	if (bc->ctr == 0 || p != 0) {
		memset(bc->buf + p, 0, (sizeof bc->buf) - p);
	}

	process_block(bc->h, bc->buf, bc->ctr, 1);
#if BLAKE2_LE
	memcpy(dst, bc->h, bc->out_len);
#else
	for (i = 0; i < 8; i ++) {
		enc32le(tmp + (i << 2), bc->h[i]);
	}
	memcpy(dst, tmp, bc->out_len);
#endif
}

/* see blake2.h */
void
blake2s(void *dst, size_t dst_len, const void *key, size_t key_len,
	const void *src, size_t src_len)
{
	blake2s_context bc;

	blake2s_init_key(&bc, dst_len, key, key_len);
	blake2s_update(&bc, src, src_len);
	blake2s_final(&bc, dst);
}

/* see blake2.h */
void
blake2s_expand(void *dst, size_t dst_len,
	const void *seed, size_t seed_len, uint64_t label)
{
	uint32_t h[8];
	uint8_t buf[64];
	size_t in_len;
	uint64_t num;

	in_len = 16 + seed_len;
	enc64le(buf, label);
	memset(buf + 8, 0, 8);
	memcpy(buf + 16, seed, seed_len);
	memset(buf + in_len, 0, (sizeof buf) - in_len);
	num = 0;
#if BLAKE2_AVX2
	if (dst_len >= 64) {
		uint8_t buf_x2[128];

		memcpy(buf_x2 +   0, buf +  0, 16);
		memcpy(buf_x2 +  16, buf +  0, 16);
		memcpy(buf_x2 +  32, buf + 16, 16);
		memcpy(buf_x2 +  48, buf + 16, 16);
		memcpy(buf_x2 +  64, buf + 32, 16);
		memcpy(buf_x2 +  80, buf + 32, 16);
		memcpy(buf_x2 +  96, buf + 48, 16);
		memcpy(buf_x2 + 112, buf + 48, 16);
		buf_x2[24] = 1;
		while (dst_len >= 64) {
			/* We compute all full pairs of blocks, but
			   only at most 2^31 pairs at a time, since
			   expand_inner_x2() cannot propagate block
			   counter carries. */
			uint64_t tnum;
			size_t tlen;

			tnum = (uint64_t)dst_len >> 6;
			if (tnum >= 0x80000000) {
				tnum = 0x80000000;
			}
			enc32le(buf_x2 + 12, (uint32_t)(num >> 32));
			enc32le(buf_x2 + 28, (uint32_t)(num >> 32));
			tlen = (size_t)tnum << 6;
			expand_inner_x2(dst, tlen, buf_x2, in_len);
			dst = (uint8_t *)dst + tlen;
			dst_len -= tlen;
			num += tnum << 1;
		}
	}
#endif
	while (dst_len > 0) {
		size_t clen;
#if !BLAKE2_LE
		uint8_t tmp[32];
		int i;
#endif

		memcpy(h, IV, sizeof h);
		h[0] ^= 0x01010000 ^ (sizeof h);
		enc64le(buf + 8, num ++);
		process_block(h, buf, in_len, 1);
		clen = dst_len < (sizeof h) ? dst_len : (sizeof h);
#if BLAKE2_LE
		memcpy(dst, h, clen);
#else
		for (i = 0; i < 8; i ++) {
			enc32le(tmp + (i << 2), h[i]);
		}
		memcpy(dst, tmp, clen);
#endif
		dst_len -= clen;
		dst = (uint8_t *)dst + clen;
	}
}
