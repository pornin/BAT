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

static const uint64_t IV[] = {
	0x6A09E667F3BCC908, 0xBB67AE8584CAA73B,
	0x3C6EF372FE94F82B, 0xA54FF53A5F1D36F1,
	0x510E527FADE682D1, 0x9B05688C2B3E6C1F,
	0x1F83D9ABFB41BD6B, 0x5BE0CD19137E2179
};

static void
process_block(uint64_t *h, const uint8_t *data, uint64_t t, int f)
{
	uint64_t v[16], m[16];
	int i;

	memcpy(v, h, 8 * sizeof(uint64_t));
	memcpy(v + 8, IV, sizeof IV);
	v[12] ^= t;
	if (f) {
		v[14] = ~v[14];
	}

#if BLAKE2_LE
	memcpy(m, data, sizeof m);
#else
	for (i = 0; i < 16; i ++) {
		m[i] = dec64le(data + (i << 3));
	}
#endif

#define ROR(x, n)   (((x) << (64 - (n))) | ((x) >> (n)))

#define G(a, b, c, d, x, y)   do { \
		v[a] += v[b] + (x); \
		v[d] = ROR(v[d] ^ v[a], 32); \
		v[c] += v[d]; \
		v[b] = ROR(v[b] ^ v[c], 24); \
		v[a] += v[b] + (y); \
		v[d] = ROR(v[d] ^ v[a], 16); \
		v[c] += v[d]; \
		v[b] = ROR(v[b] ^ v[c], 63); \
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
	ROUND( 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15);
	ROUND(14, 10,  4,  8,  9, 15, 13,  6,  1, 12,  0,  2, 11,  7,  5,  3);

#undef ROR
#undef G
#undef ROUND

	for (i = 0; i < 8; i ++) {
		h[i] ^= v[i] ^ v[i + 8];
	}
}

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
blake2b_init(blake2b_context *bc, size_t out_len)
{
	memcpy(bc->h, IV, sizeof bc->h);
	bc->h[0] ^= 0x01010000 ^ (uint64_t)out_len;
	bc->ctr = 0;
	bc->out_len = out_len;
}

/* see blake2.h */
void
blake2b_init_key(blake2b_context *bc, size_t out_len,
	const void *key, size_t key_len)
{
	blake2b_init(bc, out_len);
	if (key_len > 0) {
		bc->h[0] ^= (uint64_t)key_len << 8;
		memcpy(bc->buf, key, key_len);
		memset(bc->buf + key_len, 0, (sizeof bc->buf) - key_len);
		bc->ctr = sizeof bc->buf;
	}
}

/* see blake2.h */
void
blake2b_update(blake2b_context *bc, const void *data, size_t len)
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
blake2b_final(blake2b_context *bc, void *dst)
{
#if !BLAKE2_LE
	int i;
	uint8_t tmp[64];
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
		enc64le(tmp + (i << 3), bc->h[i]);
	}
	memcpy(dst, tmp, bc->out_len);
#endif
}

/* see blake2.h */
void
blake2b(void *dst, size_t dst_len, const void *key, size_t key_len,
	const void *src, size_t src_len)
{
	blake2b_context bc;

	blake2b_init_key(&bc, dst_len, key, key_len);
	blake2b_update(&bc, src, src_len);
	blake2b_final(&bc, dst);
}

/* see blake2.h */
void
blake2b_expand(void *dst, size_t dst_len,
	const void *seed, size_t seed_len, uint64_t label)
{
	uint64_t h[8];
	uint8_t buf[128];
	size_t in_len;
	uint64_t num;

	in_len = 16 + seed_len;
	enc64le(buf, label);
	memset(buf + 8, 0, 8);
	memcpy(buf + 16, seed, seed_len);
	memset(buf + in_len, 0, (sizeof buf) - in_len);
	num = 0;
	while (dst_len > 0) {
		size_t clen;
#if !BLAKE2_LE
		uint8_t tmp[64];
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
			enc64le(tmp + (i << 3), h[i]);
		}
		memcpy(dst, tmp, clen);
#endif
		dst_len -= clen;
		dst = (uint8_t *)dst + clen;
	}
}
