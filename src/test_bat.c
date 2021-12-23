#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "bat.h"
#include "inner.h"

static void
check_equals(const void *a1, const void *a2, size_t len, const char *msg)
{
	const uint8_t *b1, *b2;
	size_t u;

	if (memcmp(a1, a2, len) == 0) {
		return;
	}
	fprintf(stderr, "ERR: %s\n", msg);
	b1 = a1;
	b2 = a2;
	fprintf(stderr, "a1 = ");
	for (u = 0; u < len; u ++) {
		fprintf(stderr, "%02x", b1[u]);
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "a2 = ");
	for (u = 0; u < len; u ++) {
		fprintf(stderr, "%02x", b2[u]);
	}
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
}

static void
selftest_seq(uint8_t *out, size_t len, uint32_t seed)
{
	size_t i;
	uint32_t t, a, b;

	a = 0xDEAD4BAD * seed;
	b = 1;

	for (i = 0; i < len; i++) {
		t = a + b;
		a = b;
		b = t;
		out[i] = (t >> 24) & 0xFF;
	}
}

static void
test_BLAKE2s_self()
{
	/*
	 * This code is from RFC 7693 (appendix E).
	 */

	// Grand hash of hash results.
	static const uint8_t blake2s_res[32] = {
		0x6A, 0x41, 0x1F, 0x08, 0xCE, 0x25, 0xAD, 0xCD,
		0xFB, 0x02, 0xAB, 0xA6, 0x41, 0x45, 0x1C, 0xEC,
		0x53, 0xC5, 0x98, 0xB2, 0x4F, 0x4F, 0xC7, 0x87,
		0xFB, 0xDC, 0x88, 0x79, 0x7F, 0x4C, 0x1D, 0xFE
	};
	// Parameter sets.
	static const size_t b2s_md_len[4] = { 16, 20, 28, 32 };
	static const size_t b2s_in_len[6] = { 0, 3, 64, 65, 255, 1024 };

	size_t i, j, outlen, inlen;
	uint8_t in[1024], md[32], key[32];
	blake2s_context ctx;

	printf("Test BLAKE2s selftest: ");
	fflush(stdout);

	blake2s_init(&ctx, 32);

	for (i = 0; i < 4; i ++) {
		outlen = b2s_md_len[i];
		for (j = 0; j < 6; j++) {
			inlen = b2s_in_len[j];

			selftest_seq(in, inlen, inlen);
			blake2s(md, outlen, NULL, 0, in, inlen);
			blake2s_update(&ctx, md, outlen);

			selftest_seq(key, outlen, outlen);
			blake2s(md, outlen, key, outlen, in, inlen);
			blake2s_update(&ctx, md, outlen);
		}
		printf(".");
		fflush(stdout);
	}

	blake2s_final(&ctx, md);
	check_equals(md, blake2s_res, sizeof blake2s_res, "KAT");

	printf(" done.\n");
	fflush(stdout);
}

static void
test_BLAKE2s_expand(void)
{
	size_t u;

	printf("Test BLAKE2s expand: ");
	fflush(stdout);

	/* Test vector generated with python3 hashlib.blake2s()
	   implementation. */
	static const uint8_t seed[] = {
		0x4B, 0xFC, 0xB2, 0x19, 0x96, 0xAC, 0xE1, 0xE2,
		0xA1, 0xD5, 0x38, 0xC5, 0x4D, 0x10, 0x99, 0xBF,
		0x53, 0x20, 0x82, 0x62
	};

	uint64_t label = 0x4A1BE6AC1347378C;

	static const uint8_t ref[] = {
		0x78, 0x3D, 0xAA, 0x23, 0xB5, 0x2A, 0xDE, 0x32,
		0x8C, 0x44, 0xB5, 0xBF, 0x68, 0xB3, 0x8E, 0xA3,
		0x47, 0x49, 0xDB, 0x98, 0x96, 0xB4, 0xD8, 0x84,
		0xA0, 0xEB, 0xB0, 0x0B, 0x84, 0x91, 0x66, 0xBD,
		0x49, 0x56, 0x50, 0xEC, 0x3E, 0x89, 0x46, 0xF3,
		0x45, 0x26, 0xBF, 0xEA, 0x28, 0x63, 0xE3, 0x83,
		0x31, 0x64, 0xFE, 0x30, 0xE2, 0x89, 0x71, 0xFC,
		0x34, 0x0C, 0x13, 0x05, 0xBA, 0x0D, 0x51, 0x39,
		0x63, 0xD7, 0x41, 0x41, 0xCB, 0x4D, 0x74, 0xE8,
		0x3F, 0x62, 0x74, 0xA2, 0xE4, 0x12, 0xB3, 0x25,
		0x48, 0xC9, 0x3E, 0x57, 0xD3, 0x9E, 0xDD, 0xD7,
		0x7B, 0x35, 0xC4, 0xE8, 0x54, 0x2E, 0x78, 0x44,
		0xEF, 0xF0, 0x98, 0xCF, 0x82, 0x6B, 0xD0, 0x92,
		0x2E, 0xF6, 0x9E, 0xFA, 0xB3, 0x38, 0x83, 0x3B,
		0x96, 0x0C, 0xCF, 0xEA, 0xA8, 0x5E, 0xBE, 0x14,
		0x64, 0xD6, 0x35, 0xE3, 0xA8, 0x60, 0x40, 0xE5,
		0xF5, 0xEB, 0xDC, 0x55, 0xC8, 0x74, 0xEB, 0x21,
		0x19, 0x43, 0x98, 0x46, 0xCD, 0xBE, 0x22, 0x0A,
		0x0A, 0xF9, 0x07, 0xEA, 0xF0, 0xDC, 0x26, 0x80,
		0x43, 0x42, 0xA6, 0xEE, 0x4F, 0x73, 0xA3, 0x0E,
		0xB4, 0xDB, 0x10, 0x2D, 0x48, 0x8A, 0x43, 0xA9,
		0xC0, 0x8B, 0x31, 0x4F, 0x2B, 0x52, 0xBA, 0xE2,
		0x33, 0xBC, 0x32, 0xEA, 0xB7, 0xBB, 0x64, 0x2D,
		0x31, 0xDA, 0x42, 0x24, 0x7B, 0x7D, 0x34, 0x61,
		0xE3, 0x90, 0x2B, 0xA4, 0x93, 0x4A, 0x9D, 0x60,
		0x4C, 0x48, 0xA4, 0x9E, 0x27, 0x04, 0x7C, 0xE6,
		0x53, 0x12, 0x53, 0xD2, 0x8B, 0xC9, 0xCD, 0x4D,
		0x74, 0xF4, 0x96, 0x5D, 0x02, 0x37, 0xB4, 0x2D,
		0xBC, 0xAB, 0xDA, 0xEC, 0x4C, 0xE3, 0xF0, 0x57,
		0x12, 0x7F, 0xB9, 0xFD, 0xB7, 0x3A, 0xDE, 0x37,
		0xEF, 0x1B, 0x84, 0x5B, 0xFE, 0x1D, 0xEB, 0xC4,
		0x0C, 0xF9, 0xC7, 0xA7, 0xE0, 0xB6, 0xC7, 0xAB
	};

	for (u = 0; u < sizeof ref; u ++) {
		uint8_t out[1 + (sizeof ref)];

		out[u] = 0xFF;
		blake2s_expand(out, u, seed, sizeof seed, label);
		if (out[u] != 0xFF) {
			fprintf(stderr, "Output buffer overflow");
			exit(EXIT_FAILURE);
		}
		check_equals(out, ref, u, "KAT");
		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_BLAKE2b_self()
{
	/*
	 * This code is from RFC 7693 (appendix E).
	 */

	// Grand hash of hash results.
	static const uint8_t blake2b_res[32] = {
		0xC2, 0x3A, 0x78, 0x00, 0xD9, 0x81, 0x23, 0xBD,
		0x10, 0xF5, 0x06, 0xC6, 0x1E, 0x29, 0xDA, 0x56,
		0x03, 0xD7, 0x63, 0xB8, 0xBB, 0xAD, 0x2E, 0x73,
		0x7F, 0x5E, 0x76, 0x5A, 0x7B, 0xCC, 0xD4, 0x75
	};
	// Parameter sets.
	static const size_t b2b_md_len[4] = { 20, 32, 48, 64 };
	static const size_t b2b_in_len[6] = { 0, 3, 128, 129, 255, 1024 };

	size_t i, j, outlen, inlen;
	uint8_t in[1024], md[64], key[64];
	blake2b_context ctx;

	printf("Test BLAKE2b selftest: ");
	fflush(stdout);

	blake2b_init(&ctx, 32);

	for (i = 0; i < 4; i ++) {
		outlen = b2b_md_len[i];
		for (j = 0; j < 6; j++) {
			inlen = b2b_in_len[j];

			selftest_seq(in, inlen, inlen);
			blake2b(md, outlen, NULL, 0, in, inlen);
			blake2b_update(&ctx, md, outlen);

			selftest_seq(key, outlen, outlen);
			blake2b(md, outlen, key, outlen, in, inlen);
			blake2b_update(&ctx, md, outlen);
		}
		printf(".");
		fflush(stdout);
	}

	blake2b_final(&ctx, md);
	check_equals(md, blake2b_res, sizeof blake2b_res, "KAT");

	printf(" done.\n");
	fflush(stdout);
}

static void
test_BLAKE2b_expand(void)
{
	size_t u;

	printf("Test BLAKE2b expand: ");
	fflush(stdout);

	/* Test vector generated with python3 hashlib.blake2b()
	   implementation. */
	static const uint8_t seed[] = {
		0x4B, 0xFC, 0xB2, 0x19, 0x96, 0xAC, 0xE1, 0xE2,
		0xA1, 0xD5, 0x38, 0xC5, 0x4D, 0x10, 0x99, 0xBF,
		0x53, 0x20, 0x82, 0x62
	};

	uint64_t label = 0x4A1BE6AC1347378C;

	static const uint8_t ref[] = {
		0xF7, 0x50, 0xF1, 0x35, 0x88, 0x0B, 0x7F, 0xBD,
		0x1E, 0x01, 0x54, 0x42, 0x21, 0x6C, 0xAC, 0xCA,
		0x6A, 0x19, 0xF4, 0xFE, 0x76, 0xB1, 0x69, 0xF8,
		0x2B, 0xA1, 0x99, 0x14, 0x13, 0xF5, 0xB1, 0x87,
		0xD9, 0xF8, 0xA0, 0x49, 0x47, 0xF6, 0x94, 0x26,
		0x4E, 0x91, 0xF0, 0x63, 0x36, 0x56, 0x56, 0x9C,
		0x3D, 0xF2, 0xD9, 0x8D, 0x7D, 0x6D, 0x07, 0xF6,
		0x64, 0xB1, 0x25, 0x14, 0xB0, 0x80, 0xF6, 0x08,
		0x59, 0x70, 0xB0, 0xE2, 0x18, 0x2A, 0x0C, 0x9B,
		0xA6, 0x51, 0xE2, 0x73, 0xE8, 0xBF, 0x0A, 0x2F,
		0x3E, 0xD1, 0x65, 0x34, 0x95, 0x5F, 0xF1, 0x0C,
		0xB3, 0x0A, 0x45, 0xF5, 0x90, 0x71, 0x71, 0x72,
		0xCA, 0x5D, 0x58, 0x46, 0xF1, 0xDA, 0xC7, 0xE4,
		0xD4, 0x5B, 0xAE, 0x92, 0xBD, 0x6B, 0x0B, 0xA6,
		0xBF, 0xDD, 0x90, 0x24, 0x8B, 0x8B, 0xF7, 0x02,
		0x4F, 0xDB, 0x99, 0xA8, 0x42, 0x2D, 0x58, 0x51,
		0x55, 0xD5, 0xD4, 0xEA, 0x08, 0x94, 0x19, 0x99,
		0x5B, 0x25, 0xEB, 0x24, 0x48, 0x56, 0xDE, 0xEA,
		0xA7, 0x66, 0x02, 0xD8, 0x40, 0x2B, 0x3B, 0xCC,
		0x2B, 0x98, 0xA1, 0x9B, 0xEE, 0x59, 0xD2, 0x42,
		0x60, 0xF2, 0x80, 0x95, 0x4D, 0x3E, 0x93, 0xD9,
		0x17, 0x2B, 0xAF, 0x11, 0xD4, 0xE1, 0x40, 0x60,
		0x5F, 0xC9, 0x2D, 0x1D, 0xFA, 0x7F, 0x21, 0xAB,
		0x0C, 0xA2, 0xFE, 0x90, 0xD9, 0x23, 0x65, 0x52,
		0xA7, 0xE5, 0x33, 0xB6, 0xC3, 0xEA, 0xE4, 0xC0,
		0x91, 0xBA, 0x1C, 0xB5, 0x4B, 0x81, 0xAC, 0xBF,
		0xC3, 0x55, 0x82, 0xE7, 0xF9, 0x56, 0x0B, 0xD1,
		0x9F, 0x74, 0x18, 0xEB, 0x49, 0xEE, 0x55, 0x48,
		0xE6, 0x6F, 0xE6, 0x01, 0x69, 0x6A, 0x7C, 0x59,
		0x8D, 0xD0, 0x45, 0x1C, 0x14, 0x28, 0x44, 0x74,
		0x24, 0x95, 0xE0, 0xEB, 0x0A, 0x21, 0x82, 0x8D,
		0x99, 0x35, 0xC5, 0x1C, 0x68, 0x98, 0x51, 0x3A,
		0xF9, 0x7F, 0x09, 0xE7, 0xA8, 0xAB, 0x20, 0x80,
		0xCD, 0x2D, 0x46, 0x25, 0xCB, 0x7A, 0xC6, 0xC5,
		0xDC, 0xF5, 0xAC, 0x76, 0x00, 0xA0, 0xC0, 0xDA,
		0x29, 0x41, 0x5C, 0x2A, 0x0D, 0x0A, 0xE4, 0x18,
		0x73, 0x35, 0xD2, 0x8B, 0x46, 0xAA, 0x04, 0x8E,
		0x32, 0xB4, 0xA3, 0x79, 0x95, 0x0A, 0x9F, 0x4C,
		0x9F, 0x0D, 0xED, 0x67, 0xA8, 0x97, 0xEB, 0xB0,
		0xCA, 0xD9, 0xF1, 0xBB, 0x88, 0x7F, 0x14, 0xD0,
		0xD0, 0xCD, 0x7F, 0xEC, 0xAC, 0xDB, 0x7C, 0x81,
		0x3F, 0x19, 0x6C, 0x56, 0x16, 0x26, 0x4A, 0xA7,
		0xD8, 0x75, 0xC0, 0x91, 0xDA, 0x8A, 0x35, 0xDB,
		0x75, 0x34, 0x9F, 0x60, 0x57, 0x0A, 0xFD, 0xBD,
		0xBA, 0x43, 0x64, 0xB6, 0xF9, 0x63, 0x8C, 0x39,
		0x0C, 0xFF, 0x07, 0x09, 0xBB, 0xD8, 0x85, 0x19,
		0x0C, 0x2B, 0xDF, 0xF1, 0x97, 0xD7, 0xC2, 0x38,
		0x15, 0x89, 0x7A, 0x54, 0x6E, 0x6E, 0x30, 0xFC,
		0xA8, 0xD0, 0xCD, 0xC0, 0x82, 0x37, 0x0B, 0x6A,
		0x21, 0x24, 0x48, 0x85, 0x9F, 0xB3, 0xEA, 0x1B,
		0x12, 0xAF, 0x17, 0xD3, 0x20, 0x31, 0xE3, 0x35,
		0xB8, 0x78, 0xF7, 0x7B, 0x2C, 0x07, 0xAD, 0xEF,
		0x26, 0xEF, 0xCB, 0xC3, 0x59, 0x01, 0x9F, 0x73,
		0x5C, 0x88, 0xB3, 0x61, 0x6D, 0x77, 0x52, 0x30,
		0x04, 0x71, 0x28, 0xB8, 0x94, 0xF3, 0xA0, 0x30,
		0x05, 0xCD, 0x51, 0x2F, 0x90, 0x8B, 0xF1, 0x1F,
		0x52, 0xBC, 0x2B, 0x20, 0xD2, 0x52, 0xAE, 0x41,
		0x70, 0x56, 0x07, 0x84, 0x90, 0xAF, 0x3B, 0xE6,
		0xAD, 0x25, 0x11, 0x07, 0x36, 0x86, 0xFC, 0xD5,
		0xA5, 0x4A, 0xE7, 0x09, 0xBF, 0x02, 0x10, 0x82,
		0x52, 0xDB, 0x01, 0x77, 0x77, 0x2A, 0xAA, 0x3A,
		0xFD, 0x0F, 0x9E, 0x6E, 0x86, 0x0B, 0x6F, 0x77,
		0x7A, 0x5B, 0x1A, 0xD0, 0x9F, 0xFB, 0x49, 0x4B,
		0x79, 0x8D, 0x5C, 0x59, 0x9D, 0x5A, 0x0D, 0x51
	};

	for (u = 0; u < sizeof ref; u ++) {
		uint8_t out[1 + (sizeof ref)];

		out[u] = 0xFF;
		blake2b_expand(out, u, seed, sizeof seed, label);
		if (out[u] != 0xFF) {
			fprintf(stderr, "Output buffer overflow");
			exit(EXIT_FAILURE);
		}
		check_equals(out, ref, u, "KAT");
		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

/*
 * Initialize a PRNG with a given seed and extra label.
 */
static void
rand_init(prng_context *rng, const char *seed, uint64_t x)
{
	prng_init(rng, seed, strlen(seed), x);
}

/*
 * Generate a random polynomial with integer coefficients (coefficients
 * are signed and selected uniformly over num_bits).
 */
static void
rand_poly_32(prng_context *rng, int32_t *f, unsigned logn, unsigned num_bits)
{
	size_t u, n;

	n = (size_t)1 << logn;
	for (u = 0; u < n; u ++) {
		uint32_t x;

		x = (uint32_t)prng_get_u64(rng);
		f[u] = *(int32_t *)&x >> (32 - num_bits);
	}
}

static void
poly_add(int32_t *d, const int32_t *a, const int32_t *b, unsigned logn)
{
	size_t u, n;

	n = (size_t)1 << logn;
	for (u = 0; u < n; u ++) {
		d[u] = a[u] + b[u];
	}
}

static void
poly_sub(int32_t *d, const int32_t *a, const int32_t *b, unsigned logn)
{
	size_t u, n;

	n = (size_t)1 << logn;
	for (u = 0; u < n; u ++) {
		d[u] = a[u] - b[u];
	}
}

static void
poly_neg(int32_t *d, const int32_t *a, unsigned logn)
{
	size_t u, n;

	n = (size_t)1 << logn;
	for (u = 0; u < n; u ++) {
		d[u] = -a[u];
	}
}

static void
poly_mul(int32_t *d, const int32_t *a, const int32_t *b, unsigned logn)
{
	int32_t t[1024];
	size_t u, v, n;

	n = (size_t)1 << logn;
	memset(t, 0, sizeof t);
	for (u = 0; u < n; u ++) {
		for (v = 0; v < n; v ++) {
			int32_t m;

			m = a[u] * b[v];
			if ((u + v) < n) {
				t[u + v] += m;
			} else {
				t[u + v - n] -= m;
			}
		}
	}
	memcpy(d, t, n * sizeof *d);
}

static inline void
print_fnr(fnr x)
{
	fprintf(stderr, "%ld(%08lX)",
		(long)(*(int64_t *)&x.v >> 32),
		(unsigned long)(uint32_t)x.v);
}

static void
print_poly_i32(const char *name, const int32_t *a, unsigned logn)
{
	size_t u, n;

	n = (size_t)1 << logn;
	fprintf(stderr, "%s =", name);
	for (u = 0; u < n; u ++) {
		fprintf(stderr, " %ld", (long)a[u]);
	}
	fprintf(stderr, "\n");
}

static void
print_poly_fnr(const char *name, const fnr *f, unsigned logn)
{
	size_t u, n;

	n = (size_t)1 << logn;
	fprintf(stderr, "%s =", name);
	for (u = 0; u < n; u ++) {
		fprintf(stderr, " ");
		print_fnr(f[u]);
	}
	fprintf(stderr, "\n");
}

static void
check_poly_eq_round(const char *banner,
	const int32_t *a, const fnr *f, unsigned logn)
{
	size_t u, n;

	n = (size_t)1 << logn;
	for (u = 0; u < n; u ++) {
		if (fnr_round(f[u]) != a[u]) {
			break;
		}
	}
	if (u == n) {
		return;
	}
	fprintf(stderr, "ERR: %s (not equal on %zu)\n", banner, u);
	print_poly_i32("a", a, logn);
	print_poly_fnr("f", f, logn);
	fprintf(stderr, "a[%zu] = %ld\n", u, (long)a[u]);
	fprintf(stderr, "f[%zu] = ", u);
	print_fnr(f[u]);
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
}

static void
test_FFT(void)
{
	unsigned logn;

	printf("Test FFT: ");
	fflush(stdout);

	for (logn = 1; logn <= 10; logn ++) {
		prng_context rng;
		int32_t a[1024], b[1024], c[1024];
		fnr fa[1024], fb[1024], fc[1024];
		unsigned num_bits;
		size_t u, n;
		int i;

		rand_init(&rng, "test_FFT", logn);
		n = (size_t)1 << logn;

		/*
		 * If source coefficients are over k+1 bits (including
		 * sign bit), then product coefficients are at most
		 * 1+2*k+logn bits, and FFT coefficients will fit in
		 * 1+2*k+2*logn bits. We need this value to be at most 32.
		 */
		num_bits = (32 - 1 - 2 * logn) >> 1;

		for (i = 0; i < 100; i ++) {
			rand_poly_32(&rng, a, logn, 32 - logn);
			for (u = 0; u < n; u ++) {
				fa[u] = fnr_of(a[u]);
			}
			bat_FFT(fa, logn);
			bat_iFFT(fa, logn);
			check_poly_eq_round("FFT", a, fa, logn);

			rand_poly_32(&rng, a, logn, num_bits);
			rand_poly_32(&rng, b, logn, num_bits);
			for (u = 0; u < n; u ++) {
				fa[u] = fnr_of(a[u]);
				fb[u] = fnr_of(b[u]);
			}
			bat_FFT(fa, logn);
			bat_FFT(fb, logn);
			memcpy(fc, fa, n * sizeof(fnr));
			bat_poly_add(fc, fb, logn);
			bat_iFFT(fc, logn);
			poly_add(c, a, b, logn);
			check_poly_eq_round("add", c, fc, logn);

			memcpy(fc, fa, n * sizeof(fnr));
			bat_poly_sub(fc, fb, logn);
			bat_iFFT(fc, logn);
			poly_sub(c, a, b, logn);
			check_poly_eq_round("sub", c, fc, logn);

			memcpy(fc, fa, n * sizeof(fnr));
			bat_poly_neg(fc, logn);
			bat_iFFT(fc, logn);
			poly_neg(c, a, logn);
			check_poly_eq_round("neg", c, fc, logn);

			memcpy(fc, fa, n * sizeof(fnr));
			bat_poly_mul_fft(fc, fb, logn);
			bat_iFFT(fc, logn);
			poly_mul(c, a, b, logn);
			check_poly_eq_round("mul", c, fc, logn);
		}

		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
prep_tmp(void *tmp, size_t tmp_len, int i)
{
	memset(tmp, i & 0xFF, tmp_len);
}

static size_t
get_tmp_used(const void *tmp, size_t tmp_len, int i)
{
	const uint8_t *buf;
	size_t u;

	buf = tmp;
	i &= 0xFF;
	for (u = tmp_len; u > 0; u --) {
		if (buf[u - 1] != i) {
			return u;
		}
	}
	return 0;
}

static void
check_tmp_used(const char *name,
	const void *tmp, size_t tmp_len, int i, size_t max_len)
{
	size_t used_len;

	used_len = get_tmp_used(tmp, tmp_len, i);
	if (used_len > max_len) {
		fprintf(stderr, "ERR: %s: tmp usage exceeded allowance"
			" (%lu vs %lu bytes)\n",
			name,
			(unsigned long)used_len, (unsigned long)max_len);
		exit(EXIT_FAILURE);
	}
}

static void
test_kem_inner_spec(uint32_t q, unsigned logn)
{
	int i;
	union {
		uint8_t b[24 * 1024 + 8];
		uint32_t w[6 * 1024 + 2];
		uint64_t d;
	} tmp;
	int8_t f[1024], g[1024], F[1024], G[1024], G2[1024], c[1024];
	uint16_t h[1024];
	int32_t w[1024];
	int enc_fail;

	printf("[%u-%u]", (unsigned)q, 1u << logn);
	fflush(stdout);

	enc_fail = 0;

	for (i = 0; i < 100; i ++) {
		prng_context rng;
		uint8_t kg_seed[32];
		int j;

		rand_init(&rng, "kem_inner",
			((uint64_t)(q << 4 | logn) << 32) | (uint64_t)i);

		/*
		 * Generate a new key pair.
		 */
		for (;;) {
			int r;

			prep_tmp(tmp.w, sizeof tmp.w, i);
			prng_get_bytes(&rng, kg_seed, sizeof kg_seed);
			r = bat_keygen_make_fg(f, g, h, q, logn,
				kg_seed, sizeof kg_seed, tmp.w);
			check_tmp_used("bat_keygen_make_fg",
				tmp.w, sizeof tmp.w, i, 24u << logn);
			if (!r) {
				continue;
			}

			prep_tmp(tmp.w, sizeof tmp.w, i);
			r = bat_keygen_solve_FG(F, G, f, g, q, logn, tmp.w);
			check_tmp_used("bat_keygen_solve_FG",
				tmp.w, sizeof tmp.w, i, 24u << logn);
			if (!r) {
				continue;
			}

			prep_tmp(tmp.w, sizeof tmp.w, i);
			r = bat_keygen_compute_w(w, f, g, F, G, q, logn, tmp.w);
			check_tmp_used("bat_keygen_compute_w",
				tmp.w, sizeof tmp.w, i, 24u << logn);
			if (!r) {
				continue;
			}

			break;
		}

		/*
		 * Verify the key pair is behaving properly.
		 */
		prep_tmp(tmp.w, sizeof tmp.w, i);
		if (!bat_keygen_verify_FG(f, g, F, G, q, logn, tmp.w)) {
			fprintf(stderr, "bat_keygen_verify_FG() failed\n");
			exit(EXIT_FAILURE);
		}
		check_tmp_used("bat_keygen_verify_FG",
			tmp.w, sizeof tmp.w, i, 16u << logn);

		prep_tmp(tmp.w, sizeof tmp.w, i);
		if (!bat_keygen_rebuild_G(G2, f, g, F, q, logn, tmp.w)) {
			fprintf(stderr, "bat_keygen_rebuild_G() failed\n");
			exit(EXIT_FAILURE);
		}
		check_tmp_used("bat_keygen_rebuild_G",
			tmp.w, sizeof tmp.w, i, 4u << logn);
		check_equals(G, G2, 1u << logn, "rebuild G");

		/*
		 * Do some encapsulation / decapsulation.
		 */
		for (j = 0; j < 100; j ++) {
			uint8_t sbuf[128], sbuf2[128];
			int r;

			rand_init(&rng, "kem_inner_encaps",
				((uint64_t)(q << 4 | logn) << 32)
				| (uint64_t)i | ((uint64_t)j << 16));

			for (;;) {
				prng_get_bytes(&rng, sbuf, SBUF_LEN(logn));
				if (logn < 3) {
					sbuf[0] &= (1u << (1u << logn)) - 1u;
				}

				prep_tmp(tmp.w, sizeof tmp.w, i);
				switch (q) {
				case 128:
					r = bat_encrypt_128(
						c, sbuf, (const uint8_t *)h,
						logn, tmp.w);
					break;
				case 257:
					r = bat_encrypt_257(
						c, sbuf, h, logn, tmp.w);
					break;
				case 769:
					r = bat_encrypt_769(
						c, sbuf, h, logn, tmp.w);
					break;
				default:
					fprintf(stderr,
						"Unknown q: %u\n", (unsigned)q);
					exit(EXIT_FAILURE);
				}
				check_tmp_used("bat_encrypt",
					tmp.w, sizeof tmp.w, i,
					(q == 128 ? 3u : 4u) << logn);
				if (!r) {
					/*
					 * This may happen with q = 769, but
					 * not with q = 128 or q = 257.
					 */
					if (q == 769) {
						enc_fail ++;
						continue;
					}
					fprintf(stderr,
						"bat_encrypt() failed\n");
					exit(EXIT_FAILURE);
				}
				break;
			}

			prep_tmp(tmp.w, sizeof tmp.w, i);
			switch (q) {
			case 128:
				bat_decrypt_128(sbuf2,
					c, f, g, F, G, w, logn, tmp.w);
				break;
			case 257:
				bat_decrypt_257(sbuf2,
					c, f, g, F, G, w, logn, tmp.w);
				break;
			case 769:
				bat_decrypt_769(sbuf2,
					c, f, g, F, G, w, logn, tmp.w);
				break;
			default:
				fprintf(stderr, "Unknown q: %u\n", (unsigned)q);
				exit(EXIT_FAILURE);
			}
			check_tmp_used("bat_decrypt",
				tmp.w, sizeof tmp.w, i, 8u << logn);

			check_equals(sbuf, sbuf2, SBUF_LEN(logn),
				"KEM enc/dec");
		}

		printf(".");
		fflush(stdout);
	}

	printf("(%d)", enc_fail);
	fflush(stdout);
}

static void
test_kem_inner(void)
{
	unsigned logn;

	printf("Test KEM (inner):\n   ");
	fflush(stdout);
	for (logn = 1; logn <= 8; logn ++) {
		test_kem_inner_spec(128, logn);
	}
	printf("\n   ");
	fflush(stdout);
	for (logn = 1; logn <= 9; logn ++) {
		test_kem_inner_spec(257, logn);
	}
	printf("\n   ");
	fflush(stdout);
	for (logn = 1; logn <= 10; logn ++) {
		test_kem_inner_spec(769, logn);
	}
	printf("\n");
}

#define CC(x)   do { \
		int cc_err = (x); \
		if (cc_err != 0) { \
			fprintf(stderr, "%s failed with error %d\n", \
				#x, cc_err); \
			exit(EXIT_FAILURE); \
		} \
	} while (0)

static void
test_kem_128_256(void)
{
	int i;
	bat_128_256_private_key sk, sk2;
	bat_128_256_public_key pk, pk2;
	bat_128_256_ciphertext ct, ct2;
	uint8_t tmp[BAT_128_256_TMP_KEYGEN], buf[33 + 8 * 256];
	size_t len, len2;

	printf("Test KEM-128-256: ");
	fflush(stdout);

	for (i = 0; i < 100; i ++) {
		int j;

		CC(bat_128_256_keygen(&sk, tmp, sizeof tmp));

		len = bat_128_256_encode_private_key(NULL, 0, &sk, 0);
		if (len > sizeof buf) {
			fprintf(stderr, "oversized private key encoding\n");
			exit(EXIT_FAILURE);
		}
		len2 = bat_128_256_encode_private_key(buf, sizeof buf, &sk, 0);
		if (len2 != len) {
			fprintf(stderr, "private key encoding size mismatch\n");
			exit(EXIT_FAILURE);
		}

		memset(&sk2, 0, sizeof sk2);
		len2 = bat_128_256_decode_private_key(
			&sk2, buf, sizeof buf, NULL, 0);
		if (len2 != len) {
			fprintf(stderr, "private key decoding size mismatch"
				"(%zu / %zu)\n", len, len2);
			exit(EXIT_FAILURE);
		}
		check_equals(sk.seed, sk2.seed, sizeof sk.seed, "sk seed");
		check_equals(sk.f, sk2.f, sizeof sk.f, "sk f");
		check_equals(sk.g, sk2.g, sizeof sk.g, "sk g");
		check_equals(sk.F, sk2.F, sizeof sk.F, "sk F");
		check_equals(sk.G, sk2.G, sizeof sk.G, "sk G");
		check_equals(sk.w, sk2.w, sizeof sk.w, "sk w");
		check_equals(sk.h, sk2.h, sizeof sk.h, "sk h");

		len = bat_128_256_encode_private_key(NULL, 0, &sk, 1);
		if (len > sizeof buf) {
			fprintf(stderr, "oversized private key encoding"
				" (short form)\n");
			exit(EXIT_FAILURE);
		}
		len2 = bat_128_256_encode_private_key(buf, sizeof buf, &sk, 1);
		if (len2 != len) {
			fprintf(stderr, "private key encoding size mismatch"
				" (short form) (%zu vs %zu)\n", len, len2);
			exit(EXIT_FAILURE);
		}

		memset(&sk2, 0, sizeof sk2);
		len2 = bat_128_256_decode_private_key(
			&sk2, buf, sizeof buf, tmp, sizeof tmp);
		if (len2 != len) {
			fprintf(stderr, "private key decoding size mismatch"
				" (short form)\n");
			exit(EXIT_FAILURE);
		}
		check_equals(sk.seed, sk2.seed, sizeof sk.seed, "sk seed");
		check_equals(sk.f, sk2.f, sizeof sk.f, "sk f");
		check_equals(sk.g, sk2.g, sizeof sk.g, "sk g");
		check_equals(sk.F, sk2.F, sizeof sk.F, "sk F");
		check_equals(sk.G, sk2.G, sizeof sk.G, "sk G");
		check_equals(sk.w, sk2.w, sizeof sk.w, "sk w");
		check_equals(sk.h, sk2.h, sizeof sk.h, "sk h");

		bat_128_256_get_public_key(&pk, &sk);

		len = bat_128_256_encode_public_key(NULL, 0, &pk);
		if (len > sizeof buf) {
			fprintf(stderr, "oversized public key encoding\n");
			exit(EXIT_FAILURE);
		}
		len2 = bat_128_256_encode_public_key(buf, sizeof buf, &pk);
		if (len2 != len) {
			fprintf(stderr, "public key encoding size mismatch"
				" (%zu vs %zu)\n", len, len2);
			exit(EXIT_FAILURE);
		}

		memset(&pk2, 0, sizeof pk2);
		len2 = bat_128_256_decode_public_key(&pk2, buf, sizeof buf);
		if (len2 != len) {
			fprintf(stderr, "public key decoding size mismatch"
				" (%zu vs %zu)\n", len, len2);
			exit(EXIT_FAILURE);
		}
		check_equals(pk.h, pk2.h, sizeof pk.h, "pk h");

		for (j = 0; j < 100; j ++) {
			uint8_t secret[48], secret2[48];

			CC(bat_128_256_encapsulate(secret, sizeof secret,
				&ct, &pk, tmp, sizeof tmp));

			len = bat_128_256_encode_ciphertext(NULL, 0, &ct);
			if (len > sizeof buf) {
				fprintf(stderr,
					"oversized ciphertext encoding\n");
				exit(EXIT_FAILURE);
			}
			len2 = bat_128_256_encode_ciphertext(
				buf, sizeof buf, &ct);
			if (len2 != len) {
				fprintf(stderr,
					"ciphertext encoding size mismatch"
					" (%zu vs %zu)\n", len, len2);
				exit(EXIT_FAILURE);
			}

			memset(&ct2, 0, sizeof ct2);
			len2 = bat_128_256_decode_ciphertext(
				&ct2, buf, sizeof buf);
			if (len2 != len) {
				fprintf(stderr,
					"ciphertext decoding size mismatch"
					" (%zu vs %zu)\n", len, len2);
				exit(EXIT_FAILURE);
			}
			check_equals(ct.c, ct2.c, sizeof ct.c, "ct c");

			CC(bat_128_256_decapsulate(secret2, sizeof secret2,
				&ct, &sk, tmp, sizeof tmp));
			check_equals(secret, secret2, sizeof secret, "secret");
		}

		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_kem_257_512(void)
{
	int i;
	bat_257_512_private_key sk, sk2;
	bat_257_512_public_key pk, pk2;
	bat_257_512_ciphertext ct, ct2;
	uint8_t tmp[BAT_257_512_TMP_KEYGEN], buf[33 + 8 * 512];
	size_t len, len2;

	printf("Test KEM-257-512: ");
	fflush(stdout);

	for (i = 0; i < 100; i ++) {
		int j;

		CC(bat_257_512_keygen(&sk, tmp, sizeof tmp));

		len = bat_257_512_encode_private_key(NULL, 0, &sk, 0);
		if (len > sizeof buf) {
			fprintf(stderr, "oversized private key encoding\n");
			exit(EXIT_FAILURE);
		}
		len2 = bat_257_512_encode_private_key(buf, sizeof buf, &sk, 0);
		if (len2 != len) {
			fprintf(stderr, "private key encoding size mismatch\n");
			exit(EXIT_FAILURE);
		}

		memset(&sk2, 0, sizeof sk2);
		len2 = bat_257_512_decode_private_key(
			&sk2, buf, sizeof buf, NULL, 0);
		if (len2 != len) {
			fprintf(stderr, "private key decoding size mismatch"
				"(%zu / %zu)\n", len, len2);
			exit(EXIT_FAILURE);
		}
		check_equals(sk.seed, sk2.seed, sizeof sk.seed, "sk seed");
		check_equals(sk.f, sk2.f, sizeof sk.f, "sk f");
		check_equals(sk.g, sk2.g, sizeof sk.g, "sk g");
		check_equals(sk.F, sk2.F, sizeof sk.F, "sk F");
		check_equals(sk.G, sk2.G, sizeof sk.G, "sk G");
		check_equals(sk.w, sk2.w, sizeof sk.w, "sk w");
		check_equals(sk.h, sk2.h, sizeof sk.h, "sk h");

		len = bat_257_512_encode_private_key(NULL, 0, &sk, 1);
		if (len > sizeof buf) {
			fprintf(stderr, "oversized private key encoding"
				" (short form)\n");
			exit(EXIT_FAILURE);
		}
		len2 = bat_257_512_encode_private_key(buf, sizeof buf, &sk, 1);
		if (len2 != len) {
			fprintf(stderr, "private key encoding size mismatch"
				" (short form) (%zu vs %zu)\n", len, len2);
			exit(EXIT_FAILURE);
		}

		memset(&sk2, 0, sizeof sk2);
		len2 = bat_257_512_decode_private_key(
			&sk2, buf, sizeof buf, tmp, sizeof tmp);
		if (len2 != len) {
			fprintf(stderr, "private key decoding size mismatch"
				" (short form)\n");
			exit(EXIT_FAILURE);
		}
		check_equals(sk.seed, sk2.seed, sizeof sk.seed, "sk seed");
		check_equals(sk.f, sk2.f, sizeof sk.f, "sk f");
		check_equals(sk.g, sk2.g, sizeof sk.g, "sk g");
		check_equals(sk.F, sk2.F, sizeof sk.F, "sk F");
		check_equals(sk.G, sk2.G, sizeof sk.G, "sk G");
		check_equals(sk.w, sk2.w, sizeof sk.w, "sk w");
		check_equals(sk.h, sk2.h, sizeof sk.h, "sk h");

		bat_257_512_get_public_key(&pk, &sk);

		len = bat_257_512_encode_public_key(NULL, 0, &pk);
		if (len > sizeof buf) {
			fprintf(stderr, "oversized public key encoding\n");
			exit(EXIT_FAILURE);
		}
		len2 = bat_257_512_encode_public_key(buf, sizeof buf, &pk);
		if (len2 != len) {
			fprintf(stderr, "public key encoding size mismatch"
				" (%zu vs %zu)\n", len, len2);
			exit(EXIT_FAILURE);
		}

		memset(&pk2, 0, sizeof pk2);
		len2 = bat_257_512_decode_public_key(&pk2, buf, sizeof buf);
		if (len2 != len) {
			fprintf(stderr, "public key decoding size mismatch"
				" (%zu vs %zu)\n", len, len2);
			exit(EXIT_FAILURE);
		}
		check_equals(pk.h, pk2.h, sizeof pk.h, "pk h");

		for (j = 0; j < 100; j ++) {
			uint8_t secret[48], secret2[48];

			CC(bat_257_512_encapsulate(secret, sizeof secret,
				&ct, &pk, tmp, sizeof tmp));

			len = bat_257_512_encode_ciphertext(NULL, 0, &ct);
			if (len > sizeof buf) {
				fprintf(stderr,
					"oversized ciphertext encoding\n");
				exit(EXIT_FAILURE);
			}
			len2 = bat_257_512_encode_ciphertext(
				buf, sizeof buf, &ct);
			if (len2 != len) {
				fprintf(stderr,
					"ciphertext encoding size mismatch"
					" (%zu vs %zu)\n", len, len2);
				exit(EXIT_FAILURE);
			}

			memset(&ct2, 0, sizeof ct2);
			len2 = bat_257_512_decode_ciphertext(
				&ct2, buf, sizeof buf);
			if (len2 != len) {
				fprintf(stderr,
					"ciphertext decoding size mismatch"
					" (%zu vs %zu)\n", len, len2);
				exit(EXIT_FAILURE);
			}
			check_equals(ct.c, ct2.c, sizeof ct.c, "ct c");

			CC(bat_257_512_decapsulate(secret2, sizeof secret2,
				&ct, &sk, tmp, sizeof tmp));
			check_equals(secret, secret2, sizeof secret, "secret");
		}

		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

static void
test_kem_769_1024(void)
{
	int i;
	bat_769_1024_private_key sk, sk2;
	bat_769_1024_public_key pk, pk2;
	bat_769_1024_ciphertext ct, ct2;
	uint8_t tmp[BAT_769_1024_TMP_KEYGEN], buf[33 + 8 * 1024];
	size_t len, len2;

	printf("Test KEM-769-1024: ");
	fflush(stdout);

	for (i = 0; i < 100; i ++) {
		int j;

		CC(bat_769_1024_keygen(&sk, tmp, sizeof tmp));

		len = bat_769_1024_encode_private_key(NULL, 0, &sk, 0);
		if (len > sizeof buf) {
			fprintf(stderr, "oversized private key encoding\n");
			exit(EXIT_FAILURE);
		}
		len2 = bat_769_1024_encode_private_key(buf, sizeof buf, &sk, 0);
		if (len2 != len) {
			fprintf(stderr, "private key encoding size mismatch\n");
			exit(EXIT_FAILURE);
		}

		memset(&sk2, 0, sizeof sk2);
		len2 = bat_769_1024_decode_private_key(
			&sk2, buf, sizeof buf, NULL, 0);
		if (len2 != len) {
			fprintf(stderr, "private key decoding size mismatch"
				"(%zu / %zu)\n", len, len2);
			exit(EXIT_FAILURE);
		}
		check_equals(sk.seed, sk2.seed, sizeof sk.seed, "sk seed");
		check_equals(sk.f, sk2.f, sizeof sk.f, "sk f");
		check_equals(sk.g, sk2.g, sizeof sk.g, "sk g");
		check_equals(sk.F, sk2.F, sizeof sk.F, "sk F");
		check_equals(sk.G, sk2.G, sizeof sk.G, "sk G");
		check_equals(sk.w, sk2.w, sizeof sk.w, "sk w");
		check_equals(sk.h, sk2.h, sizeof sk.h, "sk h");

		len = bat_769_1024_encode_private_key(NULL, 0, &sk, 1);
		if (len > sizeof buf) {
			fprintf(stderr, "oversized private key encoding"
				" (short form)\n");
			exit(EXIT_FAILURE);
		}
		len2 = bat_769_1024_encode_private_key(buf, sizeof buf, &sk, 1);
		if (len2 != len) {
			fprintf(stderr, "private key encoding size mismatch"
				" (short form) (%zu vs %zu)\n", len, len2);
			exit(EXIT_FAILURE);
		}

		memset(&sk2, 0, sizeof sk2);
		len2 = bat_769_1024_decode_private_key(
			&sk2, buf, sizeof buf, tmp, sizeof tmp);
		if (len2 != len) {
			fprintf(stderr, "private key decoding size mismatch"
				" (short form)\n");
			exit(EXIT_FAILURE);
		}
		check_equals(sk.seed, sk2.seed, sizeof sk.seed, "sk seed");
		check_equals(sk.f, sk2.f, sizeof sk.f, "sk f");
		check_equals(sk.g, sk2.g, sizeof sk.g, "sk g");
		check_equals(sk.F, sk2.F, sizeof sk.F, "sk F");
		check_equals(sk.G, sk2.G, sizeof sk.G, "sk G");
		check_equals(sk.w, sk2.w, sizeof sk.w, "sk w");
		check_equals(sk.h, sk2.h, sizeof sk.h, "sk h");

		bat_769_1024_get_public_key(&pk, &sk);

		len = bat_769_1024_encode_public_key(NULL, 0, &pk);
		if (len > sizeof buf) {
			fprintf(stderr, "oversized public key encoding\n");
			exit(EXIT_FAILURE);
		}
		len2 = bat_769_1024_encode_public_key(buf, sizeof buf, &pk);
		if (len2 != len) {
			fprintf(stderr, "public key encoding size mismatch"
				" (%zu vs %zu)\n", len, len2);
			exit(EXIT_FAILURE);
		}

		memset(&pk2, 0, sizeof pk2);
		len2 = bat_769_1024_decode_public_key(&pk2, buf, sizeof buf);
		if (len2 != len) {
			fprintf(stderr, "public key decoding size mismatch"
				" (%zu vs %zu)\n", len, len2);
			exit(EXIT_FAILURE);
		}
		check_equals(pk.h, pk2.h, sizeof pk.h, "pk h");

		for (j = 0; j < 100; j ++) {
			uint8_t secret[48], secret2[48];

			CC(bat_769_1024_encapsulate(secret, sizeof secret,
				&ct, &pk, tmp, sizeof tmp));

			len = bat_769_1024_encode_ciphertext(NULL, 0, &ct);
			if (len > sizeof buf) {
				fprintf(stderr,
					"oversized ciphertext encoding\n");
				exit(EXIT_FAILURE);
			}
			len2 = bat_769_1024_encode_ciphertext(
				buf, sizeof buf, &ct);
			if (len2 != len) {
				fprintf(stderr,
					"ciphertext encoding size mismatch"
					" (%zu vs %zu)\n", len, len2);
				exit(EXIT_FAILURE);
			}

			memset(&ct2, 0, sizeof ct2);
			len2 = bat_769_1024_decode_ciphertext(
				&ct2, buf, sizeof buf);
			if (len2 != len) {
				fprintf(stderr,
					"ciphertext decoding size mismatch"
					" (%zu vs %zu)\n", len, len2);
				exit(EXIT_FAILURE);
			}
			check_equals(ct.c, ct2.c, sizeof ct.c, "ct c");

			CC(bat_769_1024_decapsulate(secret2, sizeof secret2,
				&ct, &sk, tmp, sizeof tmp));
			check_equals(secret, secret2, sizeof secret, "secret");
		}

		printf(".");
		fflush(stdout);
	}

	printf(" done.\n");
	fflush(stdout);
}

#if 0
/*
 * Sample code to generate key pairs and print them out in text format
 * for external analysis. Each output line contains f, g, Fd and Gd.
 */

/* defined in keygen.c */
void bat_make_Fd(int32_t *Fd, const int8_t *f, const int8_t *F,
	const int32_t *w, unsigned qp, unsigned logn, uint32_t *tmp);

static void
make_keys(uint32_t q, unsigned logn, int num)
{
	int i;
	union {
		uint8_t b[24 * 1024 + 8];
		uint32_t w[6 * 1024 + 2];
		uint64_t d;
	} tmp;
	int8_t f[1024], g[1024], F[1024], G[1024];
	uint16_t h[1024];
	int32_t w[1024], Fd[1024], Gd[1024];

	for (i = 0; i < num; i ++) {
		prng_context rng;
		uint8_t kg_seed[32];
		size_t u, n;

		rand_init(&rng, "make_keys",
			((uint64_t)(q << 4 | logn) << 32) | (uint64_t)i);

		/*
		 * Generate a new key pair.
		 */
		for (;;) {
			int r;

			prng_get_bytes(&rng, kg_seed, sizeof kg_seed);
			r = bat_keygen_make_fg(f, g, h, q, logn,
				kg_seed, sizeof kg_seed, tmp.w);
			if (!r) {
				continue;
			}

			r = bat_keygen_solve_FG(F, G, f, g, q, logn, tmp.w);
			if (!r) {
				continue;
			}

			r = bat_keygen_compute_w(w, f, g, F, G, q, logn, tmp.w);
			if (!r) {
				continue;
			}

			break;
		}
		bat_make_Fd(Fd, f, F, w, 64513, logn, tmp.w);
		bat_make_Fd(Gd, g, G, w, 64513, logn, tmp.w);

		fprintf(stderr, ".");
		fflush(stderr);

		n = 1u << logn;
		printf("[[");
		for (u = 0; u < n; u ++) {
			if (u != 0) {
				printf(", ");
			}
			printf("%d", f[u]);
		}
		printf("], [");
		for (u = 0; u < n; u ++) {
			if (u != 0) {
				printf(", ");
			}
			printf("%d", g[u]);
		}
		printf("], [");
		for (u = 0; u < n; u ++) {
			if (u != 0) {
				printf(", ");
			}
			printf("%d", Fd[u]);
		}
		printf("], [");
		for (u = 0; u < n; u ++) {
			if (u != 0) {
				printf(", ");
			}
			printf("%d", Gd[u]);
		}
		printf("]]\n");
		fflush(stdout);
	}
	fprintf(stderr, "\n");
	fflush(stderr);
}
#endif

#if 0
/*
 * Sample code to generate key pairs and export them in a custom binary
 * format for external analysis. Polynomials f, g, F and G use one byte
 * per coefficient; w uses 4 bytes per coefficient (little-endian). Keys
 * are written one after the other in the specified file.
 */

static void
make_keys_bin(const char *fname, uint32_t q, unsigned logn, int num)
{
	int i;
	union {
		uint8_t b[24 * 1024 + 8];
		uint32_t w[6 * 1024 + 2];
		uint64_t d;
	} tmp;
	int8_t f[1024], g[1024], F[1024], G[1024];
	uint16_t h[1024];
	int32_t w[1024];
	FILE *kf;

	kf = fopen(fname, "wb");
	if (kf == NULL) {
		perror("fopen");
		exit(EXIT_FAILURE);
	}
	for (i = 0; i < num; i ++) {
		prng_context rng;
		uint8_t kg_seed[32];
		size_t u, n;

		rand_init(&rng, "make_keys_bin",
			((uint64_t)(q << 4 | logn) << 32) | (uint64_t)i);

		/*
		 * Generate a new key pair.
		 */
		for (;;) {
			int r;

			prng_get_bytes(&rng, kg_seed, sizeof kg_seed);
			r = bat_keygen_make_fg(f, g, h, q, logn,
				kg_seed, sizeof kg_seed, tmp.w);
			if (!r) {
				continue;
			}

			r = bat_keygen_solve_FG(F, G, f, g, q, logn, tmp.w);
			if (!r) {
				continue;
			}

			r = bat_keygen_compute_w(w, f, g, F, G, q, logn, tmp.w);
			if (!r) {
				continue;
			}

			break;
		}

		n = 1u << logn;
		fwrite(f, 1, n, kf);
		fwrite(g, 1, n, kf);
		fwrite(F, 1, n, kf);
		fwrite(G, 1, n, kf);
		for (u = 0; u < n; u ++) {
			uint8_t tbuf[4];

			enc32le(tbuf, w[u]);
			fwrite(tbuf, 1, 4, kf);
		}

		if ((i + 1) % 100 == 0) {
			fprintf(stderr, ".");
			fflush(stderr);
		}
	}
	fprintf(stderr, "\n");
	fflush(stderr);

	fclose(kf);
}
#endif

int
main(void)
{
	test_BLAKE2s_self();
	test_BLAKE2s_expand();
	test_BLAKE2b_self();
	test_BLAKE2b_expand();

	test_FFT();
	test_kem_inner();
	test_kem_128_256();
	test_kem_257_512();
	test_kem_769_1024();

	printf("Sizes:            pub     ct   priv (short/long)\n");
	printf("BAT-128-256:     %4zu   %4zu     %3zu / %4zu\n",
		bat_128_256_encode_public_key(0, 0, 0),
		bat_128_256_encode_ciphertext(0, 0, 0),
		bat_128_256_encode_private_key(0, 0, 0, 1),
		bat_128_256_encode_private_key(0, 0, 0, 0));
	printf("BAT-257-512:     %4zu   %4zu     %3zu / %4zu\n",
		bat_257_512_encode_public_key(0, 0, 0),
		bat_257_512_encode_ciphertext(0, 0, 0),
		bat_257_512_encode_private_key(0, 0, 0, 1),
		bat_257_512_encode_private_key(0, 0, 0, 0));
	printf("BAT-769-1024:    %4zu   %4zu     %3zu / %4zu\n",
		bat_769_1024_encode_public_key(0, 0, 0),
		bat_769_1024_encode_ciphertext(0, 0, 0),
		bat_769_1024_encode_private_key(0, 0, 0, 1),
		bat_769_1024_encode_private_key(0, 0, 0, 0));

	return 0;
}
