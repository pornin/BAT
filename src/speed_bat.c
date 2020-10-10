/*
 * Speed benchmark code for BAT implementation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "bat.h"
#include "inner.h"

#ifndef DO_BENCH86
#if defined __i386__ || defined _M_IX86 || defined __x86_64__ || defined _M_X64
#define DO_BENCH86   1
#else
#define DO_BENCH86   0
#endif
#endif

#if DO_BENCH86
#include <immintrin.h>

static inline uint64_t
core_cycles(void)
{
#if defined __GNUC__ && !defined __clang__
	uint32_t hi, lo;

	_mm_lfence();
	__asm__ __volatile__ ("rdtsc" : "=d" (hi), "=a" (lo) : : );
	return ((uint64_t)hi << 32) | (uint64_t)lo;
#else
	_mm_lfence();
	return __rdtsc();
#endif
}

#endif

static void *
xmalloc(size_t len)
{
	void *buf;

	if (len == 0) {
		return NULL;
	}
	buf = malloc(len);
	if (buf == NULL) {
		fprintf(stderr, "memory allocation error\n");
		exit(EXIT_FAILURE);
	}
	return buf;
}

static void
xfree(void *buf)
{
	if (buf != NULL) {
		free(buf);
	}
}

/*
 * Benchmark function takes an opaque context and an iteration count;
 * it returns 0 on success, a negative error code on error.
 */
typedef int (*bench_fun)(void *ctx, unsigned long num);

/*
 * Returned value is the time per iteration in nanoseconds.
 * WARNING: ON x86, VALUES ARE RETURNED IN CLOCK CYCLES, NOT NANOSECONDS;
 * THRESHOLD IS IN BILLIONS OF CYCLES.
 *
 * If the benchmark function reports an error, 0.0 is returned.
 */
static double
do_bench(bench_fun bf, void *ctx, double threshold)
{
	unsigned long num;
	int r;

	/*
	 * Alsways do a few blank runs to "train" the caches and branch
	 * prediction.
	 */
	r = bf(ctx, 5);
	if (r != 0) {
		fprintf(stderr, "ERR: %d\n", r);
		return 0.0;
	}

	num = 1;
	for (;;) {
#if DO_BENCH86
		uint64_t begin, end;
#else
		clock_t begin, end;
#endif
		double tt;

#if DO_BENCH86
		begin = core_cycles();
#else
		begin = clock();
#endif
		r = bf(ctx, num);
#if DO_BENCH86
		end = core_cycles();
#else
		end = clock();
#endif
		if (r != 0) {
			fprintf(stderr, "ERR: %d\n", r);
			return 0.0;
		}
#if DO_BENCH86
		tt = (double)(end - begin) / (double)1000000000.0;
#else
		tt = (double)(end - begin) / (double)CLOCKS_PER_SEC;
#endif
		if (tt >= threshold) {
			return tt * 1000000000.0 / (double)num;
		}

		/*
		 * If the function ran for less than 0.1 seconds then
		 * we simply double the iteration number; otherwise, we
		 * use the run time to try to get a "correct" number of
		 * iterations quickly.
		 */
		if (tt < 0.1) {
			num <<= 1;
		} else {
			unsigned long num2;

			num2 = (unsigned long)((double)num
				* (threshold * 1.1) / tt);
			if (num2 <= num) {
				num2 = num + 1;
			}
			num = num2;
		}
	}
}

#define XCAT(x, y)       XCAT_(x, y)
#define XCAT_(x, y)      x ## y
#define Zn(q, n, name)   XCAT(XCAT(XCAT(bat_, q), XCAT(_, n)), XCAT(_, name))
#define Bn(q, n, name)   XCAT(XCAT(XCAT(bench_, q), XCAT(_, n)), XCAT(_, name))

#define MK_BENCH_FUNS(q, n) \
 \
typedef struct { \
	Zn(q, n, private_key) sk; \
	Zn(q, n, public_key) pk; \
	Zn(q, n, ciphertext) ct; \
	uint8_t *enc_sk; \
	size_t enc_sk_len; \
	uint8_t *enc_pk; \
	size_t enc_pk_len; \
	uint8_t *enc_ct; \
	size_t enc_ct_len; \
	uint8_t *tmp; \
	size_t tmp_len; \
	uint8_t secret[32]; \
	unsigned logn; \
	uint8_t *sbuf; \
} Bn(q, n, context); \
 \
static int \
Bn(q, n, keygen)(void *ctx, unsigned long num) \
{ \
	Bn(q, n, context) *bc; \
 \
	bc = ctx; \
	while (num -- > 0) { \
		if (Zn(q, n, keygen)(&bc->sk, bc->tmp, bc->tmp_len) != 0) { \
			return -1; \
		} \
	} \
	return 0; \
} \
 \
static int \
Bn(q, n, encode_private_key_short)(void *ctx, unsigned long num) \
{ \
	Bn(q, n, context) *bc; \
 \
	bc = ctx; \
	while (num -- > 0) { \
		if (Zn(q, n, encode_private_key)( \
			bc->enc_sk, bc->enc_sk_len, &bc->sk, 1) == 0) \
		{ \
			return -1; \
		} \
	} \
	return 0; \
} \
 \
static int \
Bn(q, n, encode_private_key_long)(void *ctx, unsigned long num) \
{ \
	Bn(q, n, context) *bc; \
 \
	bc = ctx; \
	while (num -- > 0) { \
		if (Zn(q, n, encode_private_key)( \
			bc->enc_sk, bc->enc_sk_len, &bc->sk, 0) == 0) \
		{ \
			return -1; \
		} \
	} \
	return 0; \
} \
 \
static int \
Bn(q, n, decode_private_key)(void *ctx, unsigned long num) \
{ \
	Bn(q, n, context) *bc; \
 \
	bc = ctx; \
	while (num -- > 0) { \
		if (Zn(q, n, decode_private_key)( \
			&bc->sk, bc->enc_sk, bc->enc_sk_len, \
			bc->tmp, bc->tmp_len) == 0) \
		{ \
			return -1; \
		} \
	} \
	return 0; \
} \
 \
static int \
Bn(q, n, encode_public_key)(void *ctx, unsigned long num) \
{ \
	Bn(q, n, context) *bc; \
 \
	bc = ctx; \
	while (num -- > 0) { \
		if (Zn(q, n, encode_public_key)( \
			bc->enc_pk, bc->enc_pk_len, &bc->pk) == 0) \
		{ \
			return -1; \
		} \
	} \
	return 0; \
} \
 \
static int \
Bn(q, n, decode_public_key)(void *ctx, unsigned long num) \
{ \
	Bn(q, n, context) *bc; \
 \
	bc = ctx; \
	while (num -- > 0) { \
		if (Zn(q, n, decode_public_key)( \
			&bc->pk, bc->enc_pk, bc->enc_pk_len) == 0) \
		{ \
			return -1; \
		} \
	} \
	return 0; \
} \
 \
static int \
Bn(q, n, encode_ciphertext)(void *ctx, unsigned long num) \
{ \
	Bn(q, n, context) *bc; \
 \
	bc = ctx; \
	while (num -- > 0) { \
		if (Zn(q, n, encode_ciphertext)( \
			bc->enc_ct, bc->enc_ct_len, &bc->ct) == 0) \
		{ \
			return -1; \
		} \
	} \
	return 0; \
} \
 \
static int \
Bn(q, n, decode_ciphertext)(void *ctx, unsigned long num) \
{ \
	Bn(q, n, context) *bc; \
 \
	bc = ctx; \
	while (num -- > 0) { \
		if (Zn(q, n, decode_ciphertext)( \
			&bc->ct, bc->enc_ct, bc->enc_ct_len) == 0) \
		{ \
			return -1; \
		} \
	} \
	return 0; \
} \
 \
static int \
Bn(q, n, encapsulate)(void *ctx, unsigned long num) \
{ \
	Bn(q, n, context) *bc; \
 \
	bc = ctx; \
	while (num -- > 0) { \
		if (Zn(q, n, encapsulate)( \
			bc->secret, sizeof bc->secret, \
			&bc->ct, &bc->pk, bc->tmp, bc->tmp_len) != 0) \
		{ \
			return -1; \
		} \
	} \
	return 0; \
} \
 \
static int \
Bn(q, n, decapsulate)(void *ctx, unsigned long num) \
{ \
	Bn(q, n, context) *bc; \
 \
	bc = ctx; \
	while (num -- > 0) { \
		if (Zn(q, n, decapsulate)( \
			bc->secret, sizeof bc->secret, \
			&bc->ct, &bc->sk, bc->tmp, bc->tmp_len) != 0) \
		{ \
			return -1; \
		} \
	} \
	return 0; \
} \
 \
static int \
Bn(q, n, encapsulate_nofo)(void *ctx, unsigned long num) \
{ \
	Bn(q, n, context) *bc; \
 \
	bc = ctx; \
	while (num -- > 0) { \
		if (!XCAT(bat_encapsulate_, q)(bc->ct.c, bc->sbuf, \
			bc->pk.h, bc->logn, (uint32_t *)bc->tmp)) \
		{ \
			return -1; \
		} \
	} \
	return 0; \
} \
 \
static int \
Bn(q, n, decapsulate_nofo)(void *ctx, unsigned long num) \
{ \
	Bn(q, n, context) *bc; \
 \
	bc = ctx; \
	while (num -- > 0) { \
		XCAT(bat_decapsulate_, q)(bc->sbuf, bc->ct.c, \
			bc->sk.f, bc->sk.g, bc->sk.F, bc->sk.G, \
			bc->sk.w, bc->logn, (uint32_t *)bc->tmp); \
	} \
	return 0; \
} \
 \
static void \
Bn(q, n, all)(double threshold) \
{ \
	Bn(q, n, context) bc; \
 \
	printf("q=%3u, n=%4u:", (unsigned)q, (unsigned)n); \
	fflush(stdout); \
 \
	bc.enc_sk_len = Zn(q, n, encode_private_key(0, 0, 0, 0)); \
	bc.enc_pk_len = Zn(q, n, encode_public_key(0, 0, 0)); \
	bc.enc_ct_len = Zn(q, n, encode_public_key(0, 0, 0)); \
	bc.enc_sk = xmalloc(bc.enc_sk_len); \
	bc.enc_pk = xmalloc(bc.enc_pk_len); \
	bc.enc_ct = xmalloc(bc.enc_ct_len); \
	bc.tmp_len = 24 * n + 7; \
	bc.tmp = xmalloc(bc.tmp_len); \
	for (bc.logn = 1; (1u << bc.logn) < n; bc.logn ++); \
	bc.sbuf = xmalloc(SBUF_LEN(bc.logn)); \
	if (!bat_get_seed(bc.sbuf, SBUF_LEN(bc.logn))) { \
		fprintf(stderr, "ERR: bat_get_seed() failed\n"); \
		exit(EXIT_FAILURE); \
	} \
 \
 	PRINT_BENCHS(q, n); \
 \
	xfree(bc.enc_sk); \
	xfree(bc.enc_pk); \
	xfree(bc.enc_ct); \
	xfree(bc.tmp); \
	xfree(bc.sbuf); \
}

#if DO_BENCH86
#define PRINT_BENCHS(q, n) \
	do { \
		printf(" %7.0fk", \
			do_bench(&Bn(q, n, keygen), \
			&bc, threshold) / 1000.0); \
		fflush(stdout); \
		printf(" %8.0f", \
			do_bench(&Bn(q, n, encode_private_key_short), \
			&bc, threshold)); \
		fflush(stdout); \
		printf(" %8.0f", \
			do_bench(&Bn(q, n, decode_private_key), \
			&bc, threshold)); \
		fflush(stdout); \
		printf(" %8.0f", \
			do_bench(&Bn(q, n, encode_private_key_long), \
			&bc, threshold)); \
		fflush(stdout); \
		printf(" %8.0f", \
			do_bench(&Bn(q, n, decode_private_key), \
			&bc, threshold)); \
		fflush(stdout); \
		Zn(q, n, get_public_key)(&bc.pk, &bc.sk); \
		printf(" %8.0f", \
			do_bench(&Bn(q, n, encode_public_key), \
			&bc, threshold)); \
		fflush(stdout); \
		printf(" %8.0f", \
			do_bench(&Bn(q, n, decode_public_key), \
			&bc, threshold)); \
		fflush(stdout); \
		printf(" %8.0f", \
			do_bench(&Bn(q, n, encapsulate_nofo), \
			&bc, threshold)); \
		fflush(stdout); \
		printf(" %8.0f", \
			do_bench(&Bn(q, n, encapsulate), \
			&bc, threshold)); \
		fflush(stdout); \
		printf(" %8.0f", \
			do_bench(&Bn(q, n, decapsulate_nofo), \
			&bc, threshold)); \
		fflush(stdout); \
		printf(" %8.0f", \
			do_bench(&Bn(q, n, decapsulate), \
			&bc, threshold)); \
		fflush(stdout); \
		printf(" %8.0f", \
			do_bench(&Bn(q, n, encode_ciphertext), \
			&bc, threshold)); \
		fflush(stdout); \
		printf(" %8.0f", \
			do_bench(&Bn(q, n, decode_ciphertext), \
			&bc, threshold)); \
		printf("\n"); \
		fflush(stdout); \
	} while (0)
#else
#define PRINT_BENCHS(q, n) \
	do { \
		printf(" %8.2f", \
			do_bench(&Bn(q, n, keygen), \
			&bc, threshold) / 1000000.0); \
		fflush(stdout); \
		printf(" %8.2f", \
			do_bench(&Bn(q, n, encode_private_key_short), \
			&bc, threshold) / 1000.0); \
		fflush(stdout); \
		printf(" %8.2f", \
			do_bench(&Bn(q, n, decode_private_key), \
			&bc, threshold) / 1000.0); \
		fflush(stdout); \
		printf(" %8.2f", \
			do_bench(&Bn(q, n, encode_private_key_long), \
			&bc, threshold) / 1000.0); \
		fflush(stdout); \
		printf(" %8.2f", \
			do_bench(&Bn(q, n, decode_private_key), \
			&bc, threshold) / 1000.0); \
		fflush(stdout); \
		Zn(q, n, get_public_key)(&bc.pk, &bc.sk); \
		printf(" %8.2f", \
			do_bench(&Bn(q, n, encode_public_key), \
			&bc, threshold) / 1000.0); \
		fflush(stdout); \
		printf(" %8.2f", \
			do_bench(&Bn(q, n, decode_public_key), \
			&bc, threshold) / 1000.0); \
		fflush(stdout); \
		printf(" %8.2f", \
			do_bench(&Bn(q, n, encapsulate_nofo), \
			&bc, threshold) / 1000.0); \
		fflush(stdout); \
		printf(" %8.2f", \
			do_bench(&Bn(q, n, encapsulate), \
			&bc, threshold) / 1000.0); \
		fflush(stdout); \
		printf(" %8.2f", \
			do_bench(&Bn(q, n, decapsulate_nofo), \
			&bc, threshold) / 1000.0); \
		fflush(stdout); \
		printf(" %8.2f", \
			do_bench(&Bn(q, n, decapsulate), \
			&bc, threshold) / 1000.0); \
		fflush(stdout); \
		printf(" %8.2f", \
			do_bench(&Bn(q, n, encode_ciphertext), \
			&bc, threshold) / 1000.0); \
		fflush(stdout); \
		printf(" %8.2f", \
			do_bench(&Bn(q, n, decode_ciphertext), \
			&bc, threshold) / 1000.0); \
		printf("\n"); \
		fflush(stdout); \
	} while (0)
#endif

MK_BENCH_FUNS(128, 256)
MK_BENCH_FUNS(257, 512)
MK_BENCH_FUNS(769, 1024)

int
main(int argc, char *argv[])
{
	double threshold;

	if (argc < 2) {
		threshold = 2.0;
	} else if (argc == 2) {
		threshold = atof(argv[1]);
	} else {
		threshold = -1.0;
	}
	if (threshold <= 0.0 || threshold > 60.0) {
		fprintf(stderr,
"usage: speed [ threshold ]\n"
"'threshold' is the minimum time for a bench run, in seconds (must be\n"
"positive and less than 60).\n");
		exit(EXIT_FAILURE);
	}
#if DO_BENCH86
	printf("time threshold = %.4f Gcyc\n", threshold);
#else
	printf("time threshold = %.4f s\n", threshold);
#endif
	printf("esk / dsk = encode / decode private key (s = short format, l = long format)\n");
	printf("epk / dpk = encode / decode public key\n");
	printf("ect / dct = encode / decode ciphertext\n");
	printf("ecp = encapsulate, dcp = decapsulate (nofo = without Fujisaki-Okamoto)\n");
#if DO_BENCH86
	printf("x86 PLATFORM, USING TSC; VALUES IN CLOCK CYCLES\n");
#else
	printf("keygen in milliseconds, all other times in microseconds\n");
#endif
	printf("              "
		"   keygen"
		"    esk-s"
		"    dsk-s"
		"    esk-l"
		"    dsk-l"
		"      epk"
		"      dpk"
		" ecp-nofo"
		"   ecp-fo"
		" dcp-nofo"
		"   dcp-fo"
		"      ect"
		"      dct"
		"\n");
	bench_128_256_all(threshold);
	bench_257_512_all(threshold);
	bench_769_1024_all(threshold);
	return 0;
}
