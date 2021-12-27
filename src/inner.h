#ifndef BAT_INNER_H__
#define BAT_INNER_H__

/*
 * Internal functions for BAT.
 */

/* ====================================================================== */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "blake2.h"

#if defined BAT_AVX2 && BAT_AVX2
/*
 * This implementation uses AVX2 intrinsics.
 */
#include <immintrin.h>
#ifndef BAT_LE
#define BAT_LE   1
#endif
#ifndef BAT_UNALIGNED
#define BAT_UNALIGNED   1
#endif
#if defined __GNUC__
#define TARGET_AVX2    __attribute__((target("avx2")))
#define ALIGNED_AVX2   __attribute__((aligned(32)))
#elif defined _MSC_VER && _MSC_VER
#pragma warning( disable : 4752 )
#endif
#endif

#ifndef TARGET_AVX2
#define TARGET_AVX2
#endif
#ifndef ALIGNED_AVX2
#define ALIGNED_AVX2
#endif

/*
 * Disable warning on applying unary minus on an unsigned type.
 */
#if defined _MSC_VER && _MSC_VER
#pragma warning( disable : 4146 )
#pragma warning( disable : 4244 )
#pragma warning( disable : 4267 )
#pragma warning( disable : 4334 )
#endif

/*
 * Auto-detect 64-bit architectures.
 */
#ifndef BAT_64
#if defined __x86_64__ || defined _M_X64 \
	|| defined __ia64 || defined __itanium__ || defined _M_IA64 \
	|| defined __powerpc64__ || defined __ppc64__ || defined __PPC64__ \
	|| defined __64BIT__ || defined _LP64 || defined __LP64__ \
	|| defined __sparc64__ \
	|| defined __aarch64__ || defined _M_ARM64 \
	|| defined __mips64
#define BAT_64   1
#else
#define BAT_64   0
#endif
#endif

/*
 * Auto-detect endianness and support of unaligned accesses.
 */
#if defined __i386__ || defined _M_IX86 \
	|| defined __x86_64__ || defined _M_X64 \
	|| (defined _ARCH_PWR8 \
		&& (defined __LITTLE_ENDIAN || defined __LITTLE_ENDIAN__))

#ifndef BAT_LE
#define BAT_LE   1
#endif
#ifndef BAT_UNALIGNED
#define BAT_UNALIGNED   1
#endif

#elif (defined __LITTLE_ENDIAN && __LITTLE_ENDIAN__) \
	|| (defined __BYTE_ORDER__ && defined __ORDER_LITTLE_ENDIAN__ \
		&& __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)

#ifndef BAT_LE
#define BAT_LE   1
#endif
#ifndef BAT_UNALIGNED
#define BAT_UNALIGNED   0
#endif

#else

#ifndef BAT_LE
#define BAT_LE   0
#endif
#ifndef BAT_UNALIGNED
#define BAT_UNALIGNED   0
#endif

#endif

/*
 * For seed generation:
 *
 *  - On Linux (glibc-2.25+), FreeBSD 12+ and OpenBSD, use getentropy().
 *  - On other Unix-like systems, use /dev/urandom (also a fallback for
 *    failed getentropy() calls).
 *  - On Windows, use CryptGenRandom().
 */

#ifndef BAT_RAND_GETENTROPY
#if (defined __linux && defined __GLIBC__ \
	&& (__GLIBC__ > 2 || (__GLIBC__ == 2 && __GLIBC_MINOR__ >= 25))) \
	|| (defined __FreeBSD__ && __FreeBSD__ >= 12) \
	|| defined __OpenBSD__
#define BAT_RAND_GETENTROPY   1
#else
#define BAT_RAND_GETENTROPY   0
#endif
#endif

#ifndef BAT_RAND_URANDOM
#if defined _AIX \
	|| defined __ANDROID__ \
	|| defined __FreeBSD__ \
	|| defined __NetBSD__ \
	|| defined __OpenBSD__ \
	|| defined __DragonFly__ \
	|| defined __linux__ \
	|| (defined __sun && (defined __SVR4 || defined __svr4__)) \
	|| (defined __APPLE__ && defined __MACH__)
#define BAT_RAND_URANDOM   1
#else
#define BAT_RAND_URANDOM   0
#endif
#endif

#ifndef BAT_RAND_WIN32
#if defined _WIN32 || defined _WIN64
#define BAT_RAND_WIN32   1
#else
#define BAT_RAND_WIN32   0
#endif
#endif

/*
 * Ensure all macros are defined, to avoid warnings with -Wundef.
 */
#ifndef BAT_AVX2
#define BAT_AVX2   0
#endif

/*
 * MSVC 2015 does not known the C99 keyword 'restrict'.
 */
#if defined _MSC_VER && _MSC_VER
#ifndef restrict
#define restrict   __restrict
#endif
#endif

/* ====================================================================== */
/*
 * Fixed-point numbers.
 *
 * For FFT and other computations with approximations, we use a fixed-point
 * format over 64 bits; the top 32 bits are the integral part, and the low
 * 32 bits are the fractional part.
 */

/*
 * We wrap the type into a struct in order to detect any attempt at using
 * arithmetic operators on values directly. Since all functions are inline,
 * the compiler will be able to remove the wrapper, which will then have
 * no runtime cost.
 */
typedef struct {
	uint64_t v;
} fnr;

static inline fnr
fnr_of(int32_t j)
{
	fnr x;

	x.v = (uint64_t)j << 32;
	return x;
}

static inline fnr
fnr_of_scaled32(uint64_t t)
{
	fnr x;

	x.v = t;
	return x;
}

static inline fnr
fnr_add(fnr x, fnr y)
{
	x.v += y.v;
	return x;
}

static inline fnr
fnr_sub(fnr x, fnr y)
{
	x.v -= y.v;
	return x;
}

static inline fnr
fnr_double(fnr x)
{
	x.v <<= 1;
	return x;
}

static inline fnr
fnr_neg(fnr x)
{
	x.v = (uint64_t)0 - x.v;
	return x;
}

static inline fnr
fnr_abs(fnr x)
{
	x.v -= (x.v << 1) & -(uint64_t)(x.v >> 63);
	return x;
}

static inline fnr
fnr_mul(fnr x, fnr y)
{
#if defined __GNUC__ && defined __x86_64__
	__int128 z;

	z = (__int128)*(int64_t *)&x.v * (__int128)*(int64_t *)&y.v;
	x.v = (uint64_t)(z >> 32);
	return x;
#else
	int32_t xh, yh;
	uint32_t xl, yl;
	uint64_t z0, z1, z2, z3;

	xl = (uint32_t)x.v;
	yl = (uint32_t)y.v;
	xh = (int32_t)(*(int64_t *)&x.v >> 32);
	yh = (int32_t)(*(int64_t *)&y.v >> 32);
	z0 = ((uint64_t)xl * (uint64_t)yl + 0x80000000ul) >> 32;
	z1 = (uint64_t)((int64_t)xl * (int64_t)yh);
	z2 = (uint64_t)((int64_t)yl * (int64_t)xh);
	z3 = (uint64_t)((int64_t)xh * (int64_t)yh) << 32;
	x.v = z0 + z1 + z2 + z3;
	return x;
#endif
}

static inline fnr
fnr_sqr(fnr x)
{
#if defined __GNUC__ && defined __x86_64__
	int64_t t;
	__int128 z;

	t = *(int64_t *)&x.v;
	z = (__int128)t * (__int128)t;
	x.v = (uint64_t)(z >> 32);
	return x;
#else
	int32_t xh;
	uint32_t xl;
	uint64_t z0, z1, z3;

	xl = (uint32_t)x.v;
	xh = (int32_t)(*(int64_t *)&x.v >> 32);
	z0 = ((uint64_t)xl * (uint64_t)xl + 0x80000000ul) >> 32;
	z1 = (uint64_t)((int64_t)xl * (int64_t)xh);
	z3 = (uint64_t)((int64_t)xh * (int64_t)xh) << 32;
	x.v = z0 + (z1 << 1) + z3;
	return x;
#endif
}

static inline int32_t
fnr_round(fnr x)
{
	x.v += 0x80000000ul;
	return (int32_t)(*(int64_t *)&x.v >> 32);
}

static inline fnr
fnr_div_2e(fnr x, unsigned n)
{
	int64_t v;

	v = *(int64_t *)&x.v;
	x.v = (uint64_t)((v + (((int64_t)1 << n) >> 1)) >> n);
	return x;
}

static inline fnr
fnr_mul_2e(fnr x, unsigned n)
{
	x.v <<= n;
	return x;
}

uint64_t bat_fnr_div(uint64_t x, uint64_t y);

static inline fnr
fnr_inv(fnr x)
{
	x.v = bat_fnr_div((uint64_t)1 << 32, x.v);
	return x;
}

static inline fnr
fnr_div(fnr x, fnr y)
{
	x.v = bat_fnr_div(x.v, y.v);
	return x;
}

static inline int
fnr_lt(fnr x, fnr y)
{
	return *(int64_t *)&x.v < *(int64_t *)&y.v;
}

static const fnr fnr_zero = { 0 };
static const fnr fnr_sqrt2 = { 6074001000ull };

/* ====================================================================== */
/*
 * Apply FFT on a vector.
 */
void bat_FFT(fnr *f, unsigned logn);

/*
 * Apply inverse FFT on a vector.
 */
void bat_iFFT(fnr *f, unsigned logn);

/*
 * Add polynomial b to polynomial a (works in FFT and non-FFT). The two
 * polynomial arrays must be distinct.
 */
void bat_poly_add(fnr *restrict a, const fnr *restrict b, unsigned logn);

/*
 * Subtract polynomial b from polynomial a (works in FFT and non-FFT). The two
 * polynomial arrays must be distinct.
 */
void bat_poly_sub(fnr *restrict a, const fnr *restrict b, unsigned logn);

/*
 * Negate polynomial a (works in FFT and non-FFT).
 */
void bat_poly_neg(fnr *a, unsigned logn);

/*
 * Multiply polynomial a by constant c.
 */
void bat_poly_mulconst(fnr *a, fnr c, unsigned logn);

/*
 * Multiply polynomial a by polynomial b (FFT representation only). The two
 * polynomial arrays must be distinct.
 */
void bat_poly_mul_fft(fnr *restrict a, const fnr *restrict b, unsigned logn);

/*
 * Compute the adjoint of a polynomial in FFT representation.
 */
void bat_poly_adj_fft(fnr *a, unsigned logn);

/*
 * Scale a polynomial down by a factor 2^e.
 */
void bat_poly_div_2e(fnr *a, unsigned e, unsigned logn);

/*
 * Multiply polynomial a by polynomial b (FFT representation only). The two
 * polynomial arrays must be distinct. The polynomial b must be auto-adjoint,
 * i.e. all its coefficients in FFT representation are real numbers (the
 * polynomial has half-length; the imaginary values of the coefficients,
 * assumed to be zero and located in the second half, are not accessed).
 */
void bat_poly_mul_autoadj_fft(fnr *restrict a,
	const fnr *restrict b, unsigned logn);

/*
 * Divide polynomial a by polynomial b (FFT representation only). The two
 * polynomial arrays must be distinct. The polynomial b must be auto-adjoint,
 * i.e. all its coefficients in FFT representation are real numbers (the
 * polynomial has half-length; the imaginary values of the coefficients,
 * assumed to be zero and located in the second half, are not accessed).
 */
void bat_poly_div_autoadj_fft(fnr *restrict a,
	const fnr *restrict b, unsigned logn);

/*
 * Compute (2^e)/(a*adj(a)+b*adj(b)) into d[]. Polynomials are in FFT
 * representation. d[] is a half-size polynomial because all FFT
 * coefficients are zero (they are not set by this function). Parameter e
 * can be 0.
 */
void bat_poly_invnorm_fft(fnr *restrict d,
	const fnr *restrict a, const fnr *restrict b,
	unsigned e, unsigned logn);

/* ====================================================================== */

/*
 * Max size in bits for elements of (f,g), indexed by log(N). Size includes
 * the sign bit.
 */
extern const uint8_t bat_max_fg_bits[];

/*
 * Max size in bits for elements of (F,G), indexed by log(N). Size includes
 * the sign bit.
 */
extern const uint8_t bat_max_FG_bits[];

/*
 * Max size in bits for elements of w, indexed by log(N). Size includes
 * the sign bit.
 */
extern const uint8_t bat_max_w_bits[];

/* ====================================================================== */

/*
 * Key pair generation, first step: given a seed, candidate polynomials
 * f and g are generated. The following properties are checked:
 *  - All coefficients of f and g are within the expected bounds.
 *  - Res(f, x^n+1) == 1 mod 2.
 *  - Res(g, x^n+1) == 1 mod 2.
 *  - The (f,g) vector has an acceptable norm, both in normal and in
 *    orthogonalized representations.
 *  - f is invertible modulo x^n+1 modulo q.
 * If any of these properties is not met, then a failure is reported
 * (returned value is 0) and the contents of f[] and g[] are indeterminate.
 * Otherwise, success (1) is returned.
 *
 * If h != NULL, then the public key h = g/f mod x^n+1 mod q is returned
 * in that array. Note that h is always internally computed, regardless
 * of whether h == NULL or not.
 *
 * Size of tmp[]: 6*n elements (24*n bytes).
 * tmp[] MUST be 64-bit aligned.
 *
 * The seed length MUST NOT exceed 48 bytes.
 */
int bat_keygen_make_fg(int8_t *f, int8_t *g, uint16_t *h,
	uint32_t q, unsigned logn,
	const void *seed, size_t seed_len, uint32_t *tmp);

/*
 * Given polynomials f and g, solve the NTRU equation for F and G. This
 * may fail if there is no solution, or if some intermediate value exceeds
 * an internal heuristic threshold. Returned value is 1 on success, 0
 * on failure. On failure, contents of F and G are indeterminate.
 *
 * Size of tmp[]: 6*n elements (24*n bytes).
 * tmp[] MUST be 64-bit aligned.
 */
int bat_keygen_solve_FG(int8_t *F, int8_t *G,
	const int8_t *f, const int8_t *g,
	uint32_t q, unsigned logn, uint32_t *tmp);

/*
 * Given polynomials f, g and F, rebuild the polynomial G that completes
 * the NTRU equation g*F - f*G = q. Returned value is 1 on success, 0 on
 * failure. A failure is reported if the rebuilt solution has
 * coefficients outside of the expected maximum range, or f is not
 * invertible modulo x^n+1 modulo q. This function does NOT fully verify
 * that f, g, F, G is a solution to the NTRU equation.
 *
 * Size of tmp[]: n elements (4*n bytes).
 */
int bat_keygen_rebuild_G(int8_t *G,
	const int8_t *f, const int8_t *g, const int8_t *F,
	uint32_t q, unsigned logn, uint32_t *tmp);

/*
 * Verify that the given f, g, F, G fulfill the NTRU equation g*F - f*G = q.
 * Returned value is 1 on success, 0 on error.
 *
 * This function may be called when decoding a private key of unsure
 * provenance. It is implicitly called by bat_keygen_solve_FG().
 *
 * Size of tmp[]: 4*n elements (16*n bytes).
 */
int bat_keygen_verify_FG(
	const int8_t *f, const int8_t *g, const int8_t *F, const int8_t *G,
	uint32_t q, unsigned logn, uint32_t *tmp);

/*
 * Compute the w vector. Returned value is 1 on success, 0 on error. An
 * error is reported if the w vector has coefficients that do not fit
 * in signed 16-bit integer, or if the norm of (gamma*F_d, G_d) exceeds
 * the prescribed limit.
 *
 * Size of tmp[]: 6*n elements (24*n bytes).
 * tmp[] MUST be 64-bit aligned.
 */
int bat_keygen_compute_w(int32_t *w,
	const int8_t *f, const int8_t *g, const int8_t *F, const int8_t *G,
	uint32_t q, unsigned logn, uint32_t *tmp);

/*
 * Compute the public key h = g/f. Returned value is 1 on success, 0 on
 * error. An error is reported if f is not invertible modulo X^n+1.
 * This function is for q = 128 and 1 <= logn <= 8.
 * CAUTION: for q = 128, public key is in an array of uint8_t, not uint16_t.
 *
 * Size of tmp[]: 3*n/4 elements (3*n bytes).
 */
int bat_make_public_128(uint8_t *h, const int8_t *f, const int8_t *g,
	unsigned logn, uint32_t *tmp);

/*
 * Compute the public key h = g/f. Returned value is 1 on success, 0 on
 * error. An error is reported if f is not invertible modulo X^n+1.
 * This function is for q = 257 and 1 <= logn <= 9.
 *
 * Size of tmp[]: n elements (4*n bytes).
 */
int bat_make_public_257(uint16_t *h, const int8_t *f, const int8_t *g,
	unsigned logn, uint32_t *tmp);

/*
 * Compute the public key h = g/f. Returned value is 1 on success, 0 on
 * error. An error is reported if f is not invertible modulo X^n+1.
 * This function is for q = 769 and 1 <= logn <= 10.
 *
 * Size of tmp[]: n elements (4*n bytes).
 */
int bat_make_public_769(uint16_t *h, const int8_t *f, const int8_t *g,
	unsigned logn, uint32_t *tmp);

/*
 * Given f, g and F, rebuild G, for the case q = 128. This function
 * reports a failure if (q,logn) are not supported parameters, if f is
 * not invertible modulo x^n+1 and modulo q, or if the rebuilt value G
 * has coefficients that exceed the expected maximum size.
 *
 * This function does NOT check that the returned G matches the NTRU
 * equation.
 *
 * Size of tmp[]: 3*n/4 elements (3*n bytes).
 */
int bat_rebuild_G_128(int8_t *G,
	const int8_t *f, const int8_t *g, const int8_t *F,
	unsigned logn, uint32_t *tmp);

/*
 * Given f, g and F, rebuild G, for the case q = 257. This function
 * reports a failure if (q,logn) are not supported parameters, if f is
 * not invertible modulo x^n+1 and modulo q, or if the rebuilt value G
 * has coefficients that exceed the expected maximum size.
 *
 * This function does NOT check that the returned G matches the NTRU
 * equation.
 *
 * Size of tmp[]: n elements (4*n bytes).
 */
int bat_rebuild_G_257(int8_t *G,
	const int8_t *f, const int8_t *g, const int8_t *F,
	unsigned logn, uint32_t *tmp);

/*
 * Given f, g and F, rebuild G, for the case q = 769. This function
 * reports a failure if (q,logn) are not supported parameters, if f is
 * not invertible modulo x^n+1 and modulo q, or if the rebuilt value G
 * has coefficients that exceed the expected maximum size.
 *
 * This function does NOT check that the returned G matches the NTRU
 * equation.
 *
 * Size of tmp[]: n elements (4*n bytes).
 */
int bat_rebuild_G_769(int8_t *G,
	const int8_t *f, const int8_t *g, const int8_t *F,
	unsigned logn, uint32_t *tmp);

/* ====================================================================== */

/*
 * Get the length of sbuf, for a given degree n, with n = 2^logn.
 * The logn parameter must be between 1 and 10, inclusive. Returned length
 * is in bytes, between 1 and 128, inclusive.
 */
#define SBUF_LEN(logn)   (((1 << (logn)) + 7) >> 3)

/*
 * Encrypt: given public key (in h) and secret polynomial s (in sbuf[]),
 * produce ciphertext c1 (in c).
 *
 * This function is for q = 128, with logn = 1 to 8. Ciphertext elements
 * are in the -31..+32 range. The function cannot fail, hence it always
 * returns 1.
 * CAUTION: for q = 128, public key is in an array of uint8_t, not uint16_t.
 *
 * Size of tmp[]: 3*n/4 elements (3*n bytes)
 */
uint32_t bat_encrypt_128(int8_t *c, const uint8_t *sbuf,
	const uint8_t *h, unsigned logn, uint32_t *tmp);

/*
 * Encrypt: given public key (in h) and secret polynomial s (in sbuf[]),
 * produce ciphertext c1 (in c).
 *
 * This function is for q = 257, with logn = 1 to 9. Ciphertext elements
 * are in the -64..+64 range. The function cannot fail, hence it always
 * returns 1.
 *
 * Size of tmp[]: n elements (4*n bytes).
 */
uint32_t bat_encrypt_257(int8_t *c, const uint8_t *sbuf,
	const uint16_t *h, unsigned logn, uint32_t *tmp);

/*
 * Encrypt: given public key (in h) and secret polynomial s (in sbuf[]),
 * produce ciphertext c1 (in c).
 *
 * This function is for q = 769, with logn = 1 to 10. Ciphertext elements
 * are in the -96..+96 range.
 *
 * The function may fail, if the norm of the result is too high, in which
 * case the caller should start again with a new seed (this is uncommon).
 * On failure, this function returns 0; on success, it returns 1.
 *
 * Size of tmp[]: 3*n/4 elements (3*n bytes).
 */
uint32_t bat_encrypt_769(int8_t *c, const uint8_t *sbuf,
	const uint16_t *h, unsigned logn, uint32_t *tmp);

/*
 * Decrypt: given private key (f,g,F,G,w) and ciphertext c1, extract
 * secret s. The polynomial s has length n bits (with n = 2^logn); it
 * is returned in sbuf[] (ceil(n/8) bytes; for toy versions with logn <
 * 3, the upper bits of the incomplete byte are set to zero).
 *
 * This function is for q = 128. Ciphertext elements are in the -31..+32
 * range.
 *
 * Size of tmp[]: 2*n elements (8*n bytes).
 *
 * This function never fails; for proper security, the caller must obtain
 * the message m (using the second ciphertext element c2) and check that
 * encryption of m would indeed yield exactly ciphertext c1.
 */
void bat_decrypt_128(uint8_t *sbuf, const int8_t *c,
	const int8_t *f, const int8_t *g, const int8_t *F, const int8_t *G,
	const int32_t *w, unsigned logn, uint32_t *tmp);

/*
 * Decrypt: given private key (f,g,F,G,w) and ciphertext c1, extract
 * secret s. The polynomial s has length n bits (with n = 2^logn); it
 * is returned in sbuf[] (ceil(n/8) bytes; for toy versions with logn <
 * 3, the upper bits of the incomplete byte are set to zero).
 *
 * This function is for q = 257. Ciphertext elements are in the -64..+64
 * range.
 *
 * Size of tmp[]: 2*n elements (8*n bytes).
 *
 * This function never fails; for proper security, the caller must obtain
 * the message m (using the second ciphertext element c2) and check that
 * encryption of m would indeed yield exactly ciphertext c1.
 */
void bat_decrypt_257(uint8_t *sbuf, const int8_t *c,
	const int8_t *f, const int8_t *g, const int8_t *F, const int8_t *G,
	const int32_t *w, unsigned logn, uint32_t *tmp);

/*
 * Decrypt: given private key (f,g,F,G,w) and ciphertext c1, extract
 * secret s. The polynomial s has length n bits (with n = 2^logn); it
 * is returned in sbuf[] (ceil(n/8) bytes; for toy versions with logn <
 * 3, the upper bits of the incomplete byte are set to zero).
 *
 * This function is for q = 769. Ciphertext elements are in the -96..+96
 * range.
 *
 * Size of tmp[]: 2*n elements (8*n bytes).
 *
 * This function never fails; for proper security, the caller must obtain
 * the message m (using the second ciphertext element c2) and check that
 * encryption of m would indeed yield exactly ciphertext c1.
 */
void bat_decrypt_769(uint8_t *sbuf, const int8_t *c,
	const int8_t *f, const int8_t *g, const int8_t *F, const int8_t *G,
	const int32_t *w, unsigned logn, uint32_t *tmp);

/*
 * Second phase of decapsulation, performed modulo 769.
 * Given c', c'', f, F and w, this function computes:
 *    Fd = q'*F - f*w
 *    q*q'*Q*s' = Fd*c' - f*c''
 *
 * On input, cp[] and cs[] must contain c' and c'', respectively, in
 * Montgomery representation modulo 769. On output, polynomial q*q'*Q*s'
 * is returned in cp[], in Montgomery representation modulo 769 (since
 * coefficients of s' can have only a few specific values, this is enough
 * to recover s'). cs[] is consumed. tmp[] must have room for 4*n bytes
 * (n 32-bit elements).
 *
 * Size of tmp[]: n elements (4*n bytes).
 */
void bat_finish_decapsulate_769(uint16_t *cp, uint16_t *cs,
	const int8_t *f, const int8_t *F, const int32_t *w, unsigned logn,
	uint32_t *tmp);

/*
 * Second phase of decapsulation, performed modulo 257.
 * Given c', c'', f, F and w, this function computes:
 *    Fd = q'*F - f*w
 *    q*q'*Q*s' = Fd*c' - f*c''
 *
 * On input, cp[] and cs[] must contain c' and c'', respectively, in
 * Montgomery representation modulo 257. On output, polynomial q*q'*Q*s'
 * is returned in cp[], in Montgomery representation modulo 257 (since
 * coefficients of s' can have only a few specific values, this is enough
 * to recover s'). cs[] is consumed. tmp[] must have room for 4*n bytes
 * (n 32-bit elements).
 *
 * Size of tmp[]: n elements (4*n bytes).
 */
void bat_finish_decapsulate_257(uint16_t *cp, uint16_t *cs,
	const int8_t *f, const int8_t *F, const int32_t *w, unsigned logn,
	uint32_t *tmp);

/*
 * Explicit reduction and conversion to Montgomery representation modulo
 * 257. This works for inputs x in range 0..4278190336.
 */
static inline uint32_t
m257_tomonty(uint32_t x)
{
	x *= 16711935;
	x = (x >> 16) * 257;
	return (x >> 16) + 1;
}

/*
 * Explicit reduction and conversion to Montgomery representation modulo
 * 769. This works for inputs x in range 0..4244636416.
 */
static inline uint32_t
m769_tomonty(uint32_t x)
{
	x *= 452395775;
	x = (x >> 16) * 769;
	x = (x >> 16) + 1;
	x *= 2016233021;
	x = (x >> 16) * 769;
	return (x >> 16) + 1;
}

/* ====================================================================== */
/*
 * Computations on polynomials modulo q' = 64513.
 */

/*
 * Compute d = -a*b mod X^n+1 mod q'
 * Coefficients of source values are plain integers (for value b, they must
 * be in the -503109..+503109 range). Coefficients of output values are
 * normalized in -32256..+32256.
 *
 * Array d[] may overlap, partially or totally, with a[]; however, it
 * MUST NOT overlap with b[].
 *
 * Size of tmp[]: n/2 elements (2*n bytes).
 */
void bat_polyqp_mulneg(int16_t *d, const int16_t *a, const int32_t *b,
	unsigned logn, uint32_t *tmp);

/* ====================================================================== */
/*
 * Encoding/decoding functions.
 */

#if BAT_LE && BAT_UNALIGNED

static inline unsigned
dec16le(const void *src)
{
	return *(const uint16_t *)src;
}

static inline void
enc16le(void *dst, unsigned x)
{
	*(uint16_t *)dst = x;
}

static inline uint32_t
dec32le(const void *src)
{
	return *(const uint32_t *)src;
}

static inline void
enc32le(void *dst, uint32_t x)
{
	*(uint32_t *)dst = x;
}

static inline uint64_t
dec64le(const void *src)
{
	return *(const uint64_t *)src;
}

static inline void
enc64le(void *dst, uint64_t x)
{
	*(uint64_t *)dst = x;
}

#else

static inline unsigned
dec16le(const void *src)
{
	const uint8_t *buf;

	buf = src;
	return (unsigned)buf[0]
		| ((unsigned)buf[1] << 8);
}

static inline void
enc16le(void *dst, unsigned x)
{
	uint8_t *buf;

	buf = dst;
	buf[0] = (uint8_t)x;
	buf[1] = (uint8_t)(x >> 8);
}

static inline uint32_t
dec32le(const void *src)
{
	const uint8_t *buf;

	buf = src;
	return (uint32_t)buf[0]
		| ((uint32_t)buf[1] << 8)
		| ((uint32_t)buf[2] << 16)
		| ((uint32_t)buf[3] << 24);
}

static inline void
enc32le(void *dst, uint32_t x)
{
	uint8_t *buf;

	buf = dst;
	buf[0] = (uint8_t)x;
	buf[1] = (uint8_t)(x >> 8);
	buf[2] = (uint8_t)(x >> 16);
	buf[3] = (uint8_t)(x >> 24);
}

static inline uint64_t
dec64le(const void *src)
{
	const uint8_t *buf;

	buf = src;
	return (uint64_t)buf[0]
		| ((uint64_t)buf[1] << 8)
		| ((uint64_t)buf[2] << 16)
		| ((uint64_t)buf[3] << 24)
		| ((uint64_t)buf[4] << 32)
		| ((uint64_t)buf[5] << 40)
		| ((uint64_t)buf[6] << 48)
		| ((uint64_t)buf[7] << 56);
}

static inline void
enc64le(void *dst, uint64_t x)
{
	uint8_t *buf;

	buf = dst;
	buf[0] = (uint64_t)x;
	buf[1] = (uint64_t)(x >> 8);
	buf[2] = (uint64_t)(x >> 16);
	buf[3] = (uint64_t)(x >> 24);
	buf[4] = (uint64_t)(x >> 32);
	buf[5] = (uint64_t)(x >> 40);
	buf[6] = (uint64_t)(x >> 48);
	buf[7] = (uint64_t)(x >> 56);
}

#endif

static inline uint32_t
dec24le(const void *src)
{
	const uint8_t *buf;

	buf = src;
	return (uint32_t)buf[0]
		| ((uint32_t)buf[1] << 8)
		| ((uint32_t)buf[2] << 16);
}

static inline void
enc24le(void *dst, uint32_t x)
{
	uint8_t *buf;

	buf = dst;
	buf[0] = (uint8_t)x;
	buf[1] = (uint8_t)(x >> 8);
	buf[2] = (uint8_t)(x >> 16);
}

/*
 * bat_trim_i32_encode() and bat_trim_i32_decode() encode and decode
 * polynomials with signed coefficients (int32_t), using the specified
 * number of bits for each coefficient. The number of bits includes the
 * sign bit. Each coefficient x must be such that |x| < 2^(bits-1) (the
 * value -2^(bits-1), though conceptually encodable with two's
 * complement representation, is forbidden).
 *
 * bat_trim_i8_encode() and bat_trim_i8_decode() do the same work for
 * polynomials whose coefficients are held in slots of type int8_t.
 *
 * Encoding API:
 *
 *   Output buffer (out[]) has max length max_out_len (in bytes). If
 *   that length is not large enough, then no encoding occurs and the
 *   function returns 0; otherwise, the function returns the number of
 *   bytes which have been written into out[]. If out == NULL, then
 *   max_out_len is ignored, and no output is produced, but the function
 *   returns how many bytes it would produce.
 *
 *   Encoding functions assume that the input is valid (all values in
 *   the encodable range).
 *
 * Decoding API:
 *
 *   Input buffer (in[]) has maximum length max_in_len (in bytes). If
 *   the input length is not enough for the expected polynomial, then
 *   no decoding occurs and the function returns 0. Otherwise, the values
 *   are decoded and the number of processed input bytes is returned.
 *
 *   If the input is invalid in some way (a decoded coefficient has
 *   value -2^(bits-1), or some of the ignored bits in the last byte
 *   are non-zero), then the function fails and returns 0; the contents
 *   of the output array are then indeterminate.
 *
 * Both encoding and decoding are constant-time with regards to the
 * values and bits.
 */

size_t bat_trim_i32_encode(void *out, size_t max_out_len,
	const int32_t *x, unsigned logn, unsigned bits);
size_t bat_trim_i32_decode(int32_t *x, unsigned logn, unsigned bits,
	const void *in, size_t max_in_len);
size_t bat_trim_i8_encode(void *out, size_t max_out_len,
	const int8_t *x, unsigned logn, unsigned bits);
size_t bat_trim_i8_decode(int8_t *x, unsigned logn, unsigned bits,
	const void *in, size_t max_in_len);

/*
 * Encode a polynomial with coefficients modulo 128. This is used for
 * public keys with q = 128.
 *
 * If out == NULL, then max_out_len is ignored and the function returns
 * the size of the output it could produce (in bytes).
 * If out != NULL, then max_out_len is compared with the expected output
 * size. If max_out_len is lower, then no output is produced, and the
 * function returns 0; otherwise, the output is produced and its length
 * (in bytes) is returned.
 */
size_t bat_encode_128(void *out, size_t max_out_len,
	const uint8_t *x, unsigned logn);

/*
 * Decode a polynomial with coefficients modulo 128. This is used for
 * public keys with q = 128.
 *
 * Input buffer (in[]) has maximum length max_in_len (in bytes). If
 * the input length is not enough for the expected polynomial, then
 * no decoding occurs and the function returns 0. Otherwise, the values
 * are decoded and the number of processed input bytes is returned.
 *
 * If the input is invalid in some way (a decoded coefficient is out of
 * the expected range, or some ignored bit is non-zero), then this the
 * function fails and returns 0; the contents of the output array are
 * then indeterminate.
 *
 * Decoding is constant-time as long as no failure occurs.
 */
size_t bat_decode_128(uint8_t *x, unsigned logn,
	const void *in, size_t max_in_len);

/*
 * Encode a ciphertext polynomial, for q = 128; coefficients are in -31..+32.
 *
 * If out == NULL, then max_out_len is ignored and the function returns
 * the size of the output it could produce (in bytes).
 * If out != NULL, then max_out_len is compared with the expected output
 * size. If max_out_len is lower, then no output is produced, and the
 * function returns 0; otherwise, the output is produced and its length
 * (in bytes) is returned.
 */
size_t bat_encode_ciphertext_128(void *out, size_t max_out_len,
	const int8_t *c, unsigned logn);
/*
 * Decode a ciphertext polynomial, for q = 128; coefficients are in -31..+32.
 *
 * Input buffer (in[]) has maximum length max_in_len (in bytes). If
 * the input length is not enough for the expected polynomial, then
 * no decoding occurs and the function returns 0. Otherwise, the values
 * are decoded and the number of processed input bytes is returned.
 *
 * If the input is invalid in some way (a decoded coefficient is out of
 * the expected range, or some ignored bit is non-zero), then this the
 * function fails and returns 0; the contents of the output array are
 * then indeterminate.
 *
 * Decoding is constant-time with regard to the coefficient values.
 */
size_t bat_decode_ciphertext_128(int8_t *c, unsigned logn,
	const void *in, size_t max_in_len);

/*
 * Encode a polynomial with coefficients modulo 257. This is used for
 * public keys with q = 257.
 *
 * If out == NULL, then max_out_len is ignored and the function returns
 * the size of the output it could produce (in bytes).
 * If out != NULL, then max_out_len is compared with the expected output
 * size. If max_out_len is lower, then no output is produced, and the
 * function returns 0; otherwise, the output is produced and its length
 * (in bytes) is returned.
 */
size_t bat_encode_257(void *out, size_t max_out_len,
	const uint16_t *x, unsigned logn);

/*
 * Decode a polynomial with coefficients modulo 257. This is used for
 * public keys with q = 257.
 *
 * Input buffer (in[]) has maximum length max_in_len (in bytes). If
 * the input length is not enough for the expected polynomial, then
 * no decoding occurs and the function returns 0. Otherwise, the values
 * are decoded and the number of processed input bytes is returned.
 *
 * If the input is invalid in some way (a decoded coefficient is out of
 * the expected range, or some ignored bit is non-zero), then this the
 * function fails and returns 0; the contents of the output array are
 * then indeterminate.
 *
 * Decoding is constant-time as long as no failure occurs.
 */
size_t bat_decode_257(uint16_t *x, unsigned logn,
	const void *in, size_t max_in_len);

/*
 * Encode a ciphertext polynomial, for q = 257; coefficients are in -64..+64.
 *
 * If out == NULL, then max_out_len is ignored and the function returns
 * the size of the output it could produce (in bytes).
 * If out != NULL, then max_out_len is compared with the expected output
 * size. If max_out_len is lower, then no output is produced, and the
 * function returns 0; otherwise, the output is produced and its length
 * (in bytes) is returned.
 */
size_t bat_encode_ciphertext_257(void *out, size_t max_out_len,
	const int8_t *c, unsigned logn);
/*
 * Decode a ciphertext polynomial, for q = 257; coefficients are in -64..+64.
 *
 * Input buffer (in[]) has maximum length max_in_len (in bytes). If
 * the input length is not enough for the expected polynomial, then
 * no decoding occurs and the function returns 0. Otherwise, the values
 * are decoded and the number of processed input bytes is returned.
 *
 * If the input is invalid in some way (a decoded coefficient is out of
 * the expected range, or some ignored bit is non-zero), then this the
 * function fails and returns 0; the contents of the output array are
 * then indeterminate.
 *
 * Decoding is constant-time with regard to the coefficient values.
 */
size_t bat_decode_ciphertext_257(int8_t *c, unsigned logn,
	const void *in, size_t max_in_len);

/*
 * Encode a polynomial with coefficients modulo 769. This is used for
 * public keys with q = 769.
 *
 * If out == NULL, then max_out_len is ignored and the function returns
 * the size of the output it could produce (in bytes).
 * If out != NULL, then max_out_len is compared with the expected output
 * size. If max_out_len is lower, then no output is produced, and the
 * function returns 0; otherwise, the output is produced and its length
 * (in bytes) is returned.
 */
size_t bat_encode_769(void *out, size_t max_out_len,
	const uint16_t *x, unsigned logn);

/*
 * Decode a polynomial with coefficients modulo 769. This is used for
 * public keys with q = 769.
 *
 * Input buffer (in[]) has maximum length max_in_len (in bytes). If
 * the input length is not enough for the expected polynomial, then
 * no decoding occurs and the function returns 0. Otherwise, the values
 * are decoded and the number of processed input bytes is returned.
 *
 * If the input is invalid in some way (a decoded coefficient is out of
 * the expected range, or some ignored bit is non-zero), then this the
 * function fails and returns 0; the contents of the output array are
 * then indeterminate.
 *
 * Decoding is constant-time with regard to the coefficient values.
 */
size_t bat_decode_769(uint16_t *x, unsigned logn,
	const void *in, size_t max_in_len);

/*
 * Encode a ciphertext polynomial, for q = 769; coefficients are in -96..+96.
 *
 * If out == NULL, then max_out_len is ignored and the function returns
 * the size of the output it could produce (in bytes).
 * If out != NULL, then max_out_len is compared with the expected output
 * size. If max_out_len is lower, then no output is produced, and the
 * function returns 0; otherwise, the output is produced and its length
 * (in bytes) is returned.
 */
size_t bat_encode_ciphertext_769(void *out, size_t max_out_len,
	const int8_t *c, unsigned logn);
/*
 * Decode a ciphertext polynomial, for q = 769; coefficients are in -96..+96.
 *
 * Input buffer (in[]) has maximum length max_in_len (in bytes). If
 * the input length is not enough for the expected polynomial, then
 * no decoding occurs and the function returns 0. Otherwise, the values
 * are decoded and the number of processed input bytes is returned.
 *
 * If the input is invalid in some way (a decoded coefficient is out of
 * the expected range, or some ignored bit is non-zero), then this the
 * function fails and returns 0; the contents of the output array are
 * then indeterminate.
 *
 * Decoding is constant-time with regard to the coefficient values.
 */
size_t bat_decode_ciphertext_769(int8_t *c, unsigned logn,
	const void *in, size_t max_in_len);

/* ====================================================================== */

/*
 * Obtain a random seed from the system RNG. Maximum allowed seed length
 * is 2048 bits (256 bytes).
 *
 * Returned value is 1 on success, 0 on error.
 */
int bat_get_seed(void *seed, size_t len);

/*
 * Custom PRNG that outputs 64-bit integers. It is based on BLAKE2s.
 */
typedef struct {
	uint8_t buf[128];
	uint8_t key[32];
	uint64_t ctr;
	size_t ptr;
} prng_context;

/*
 * Initialize the PRNG from the provided seed and an extra 64-bit integer.
 * The seed length MUST NOT exceed 48 bytes.
 */
static inline void
prng_init(prng_context *p, const void *seed, size_t seed_len, uint64_t label)
{
	blake2s_expand(p->key, sizeof p->key, seed, seed_len, label);
	p->ctr = 0;
	p->ptr = sizeof p->buf;
}

/*
 * Get a 64-bit integer out of a PRNG.
 */
static inline uint64_t
prng_get_u64(prng_context *p)
{
	uint64_t x;

	if (p->ptr == sizeof p->buf) {
		blake2s_expand(p->buf, sizeof p->buf,
			p->key, sizeof p->key, p->ctr ++);
		p->ptr = 0;
	}
	x = dec64le(p->buf + p->ptr);
	p->ptr += 8;
	return x;
}

/*
 * Get arbitrary bytes out of a PRNG.
 */
static inline void
prng_get_bytes(prng_context *p, void *dst, size_t len)
{
	blake2s_expand(dst, len, p->key, sizeof p->key, p->ctr ++);
}

/* ====================================================================== */

#endif
