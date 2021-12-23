#ifndef BLAKE2_H__
#define BLAKE2_H__

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	uint8_t buf[64];
	uint32_t h[8];
	uint64_t ctr;
	size_t out_len;
} blake2s_context;

void blake2s_init(blake2s_context *bc, size_t out_len);

void blake2s_init_key(blake2s_context *bc, size_t out_len,
	const void *key, size_t key_len);

void blake2s_update(blake2s_context *bc, const void *data, size_t len);

void blake2s_final(blake2s_context *bc, void *dst);

void blake2s(void *dst, size_t dst_len, const void *key, size_t key_len,
	const void *src, size_t src_len);

/*
 * Use BLAKE2s as a PRNG: for a given seed, compute the concatenation:
 *   H(label || 0 || seed) || H(label || 1 || seed) || ...
 * with:
 *   H = BLAKE2s with a 32-byte output
 *   seed = provided seed (length MUST be at most 48 bytes)
 *   label = provided value (64-bit, little-endian)
 *   0, 1,... = block counter (64-bit, little-endian)
 * The concatenation output is truncated to dst_len and written in dst[].
 * The seed and dst buffers may overlap arbitrarily.
 */
void blake2s_expand(void *dst, size_t dst_len,
	const void *seed, size_t seed_len, uint64_t label);

typedef struct {
	uint8_t buf[128];
	uint64_t h[8];
	uint64_t ctr;
	size_t out_len;
} blake2b_context;

void blake2b_init(blake2b_context *bc, size_t out_len);

void blake2b_init_key(blake2b_context *bc, size_t out_len,
	const void *key, size_t key_len);

void blake2b_update(blake2b_context *bc, const void *data, size_t len);

void blake2b_final(blake2b_context *bc, void *dst);

void blake2b(void *dst, size_t dst_len, const void *key, size_t key_len,
	const void *src, size_t src_len);

/*
 * Use BLAKE2b as a PRNG: for a given seed, compute the concatenation:
 *   H(label || 0 || seed) || H(label || 1 || seed) || ...
 * with:
 *   H = BLAKE2b with a 64-byte output
 *   seed = provided seed (length MUST be at most 112 bytes)
 *   label = provided value (64-bit, little-endian)
 *   0, 1,... = block counter (64-bit, little-endian)
 * The concatenation output is truncated to dst_len and written in dst[].
 * The seed and dst buffers may overlap arbitrarily.
 */
void blake2b_expand(void *dst, size_t dst_len,
	const void *seed, size_t seed_len, uint64_t label);

#ifdef __cplusplus
}
#endif

#endif
