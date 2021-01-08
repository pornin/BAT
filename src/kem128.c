#include "inner.h"

/*
 * We use computations modulo 256 (usually implicitly through the use
 * of uint8_t as a storage type). When a value is needed modulo 128,
 * we apply a mask explicitly.
 *
 * We can use int8_t* and uint8_t* interchangeably, since the C standard
 * guarantees two's-complement and compatibility of formats for
 * exact-width types.
 *
 * Note that invertibility modulo 256 is equivalent to invertibility
 * modulo 128, since this boils down to the parity of the value at the
 * deepest recursion level (see mq_poly_inv_inner() for details). In BAT,
 * the keygen makes f with odd parity only (this is required for the
 * NTRU solving algorithm), thus f is always invertible modulo 128 (and
 * 256).
 *
 * Functions use a logarithm "stride" for access: successive elements of
 * polynomial a[] are a[0], a[1 << ls], a[2 << ls], ...
 */

static void
mq_poly_add_inner(uint8_t *d,
	const uint8_t *a, const uint8_t *b, int ls, unsigned logn)
{
	size_t u, n;

	n = (size_t)1 << logn;
	for (u = 0; u < n; u ++) {
		d[u << ls] = a[u << ls] + b[u << ls];
	}
}

static void
mq_poly_sub_inner(uint8_t *d,
	const uint8_t *a, const uint8_t *b, int ls, unsigned logn)
{
	size_t u, n;

	n = (size_t)1 << logn;
	for (u = 0; u < n; u ++) {
		d[u << ls] = a[u << ls] - b[u << ls];
	}
}

static void
mq_poly_neg_inner(uint8_t *d, const uint8_t *a, int ls, unsigned logn)
{
	size_t u, n;

	n = (size_t)1 << logn;
	for (u = 0; u < n; u ++) {
		d[u << ls] = -a[u << ls];
	}
}

/*
 * d <- a + X*b
 */
static void
mq_poly_add_mulX_inner(uint8_t *d,
	const uint8_t *a, const uint8_t *b, int ls, unsigned logn)
{
	/*
	 * We must take care to perform the loop in a way that does not
	 * break when d == b.
	 */
	size_t u, n;
	int t;

	n = (size_t)1 << logn;
	t = -b[(n - 1) << ls];
	for (u = 0; u < n; u ++) {
		int tn;

		tn = b[u << ls];
		d[u << ls] = a[u << ls] + t;
		t = tn;
	}
}

/*
 * d <- a - X*b
 */
static void
mq_poly_sub_mulX_inner(uint8_t *d,
	const uint8_t *a, const uint8_t *b, int ls, unsigned logn)
{
	/*
	 * We must take care to perform the loop in a way that does not
	 * break when d == b.
	 */
	size_t u, n;
	int t;

	n = (size_t)1 << logn;
	t = -b[(n - 1) << ls];
	for (u = 0; u < n; u ++) {
		int tn;

		tn = b[u << ls];
		d[u << ls] = a[u << ls] - t;
		t = tn;
	}
}

/*
 * For multiplications, we use Karatsuba, with even/odd splits:
 *
 *   a = a_e(X^2) + X*a_o(X^2)
 *   b = b_e(X^2) + X*b_o(X^2)
 *   a*b = (a_e*b_e + X*a_o*b_o)(X^2) + X*(a_e*b_o + a_o*b_e)(X^2)
 *   (a_e*b_o + a_o*b_e) = (a_e + a_o)*(b_e + b_o) - a_e*b_e - a_o*b_o
 *
 * Size of tmp[]: 2*n bytes (with n = top-level degree).
 */
static void
mq_poly_mul_inner(uint8_t *d, const uint8_t *a, const uint8_t *b,
	int ls, unsigned logn, uint8_t *tmp)
{
	uint8_t *t1, *t2;

	switch (logn) {
		unsigned a0, a1, a2, a3;
		unsigned b0, b1, b2, b3;

	case 1:
		a0 = a[0 << ls];
		a1 = a[1 << ls];
		b0 = b[0 << ls];
		b1 = b[1 << ls];
		d[0 << ls] = a0 * b0 - a1 * b1;
		d[1 << ls] = a0 * b1 + a1 * b0;
		return;
	case 2:
		a0 = a[0 << ls];
		a1 = a[1 << ls];
		a2 = a[2 << ls];
		a3 = a[3 << ls];
		b0 = b[0 << ls];
		b1 = b[1 << ls];
		b2 = b[2 << ls];
		b3 = b[3 << ls];
		d[0 << ls] = a0 * b0 - a1 * b3 - a2 * b2 - a3 * b1;
		d[1 << ls] = a0 * b1 + a1 * b0 - a2 * b3 - a3 * b2;
		d[2 << ls] = a0 * b2 + a1 * b1 + a2 * b0 - a3 * b3;
		d[3 << ls] = a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0;
		return;
	default:
		break;
	}

	/*
	 * a_e is a[], starting at 0, with stride ls + 1
	 * a_o is a[], starting at 1 << ls, with stride ls + 1
	 *
	 * We need two temporaries t1 and t2, and we will use only the
	 * elements with stride ls + 1 in them. Thus, in the received
	 * tmp[], we use odd-indexed elements for our temporaries,
	 * leaving the even-indexed elements free for the deeper
	 * recursion levels.
	 */
	t1 = tmp + (1 << ls);
	t2 = t1 + ((size_t)1 << (logn + ls));

	/*
	 * t1 <- a_e + a_o
	 * t2 <- b_e + b_o
	 */
	mq_poly_add_inner(t1, a, a + (1 << ls), ls + 1, logn - 1);
	mq_poly_add_inner(t2, b, b + (1 << ls), ls + 1, logn - 1);

	/*
	 * t1 <- (a_e + a_o)*(b_e + b_o)
	 */
	mq_poly_mul_inner(t1, t1, t2, ls + 1, logn - 1, tmp);

	/*
	 * t2 <- a_o * b_o
	 * d_e <- a_e * b_e
	 * We don't need a[] and b[] afterwards, which is why we can
	 * write into d_e (which may overlap either or both).
	 */
	mq_poly_mul_inner(t2,
		a + (1 << ls), b + (1 << ls), ls + 1, logn - 1, tmp);
	mq_poly_mul_inner(d,
		a, b, ls + 1, logn - 1, tmp);

	/*
	 * d_o <- t1 - t2 - d_e = a_e*b_o + a_o*b_e
	 */
	mq_poly_sub_inner(t1, t1, t2, ls + 1, logn - 1);
	mq_poly_sub_inner(d + (1 << ls), t1, d, ls + 1, logn - 1);

	/*
	 * d_e <- d_e + X*t2 = a_e*b_e + X*a_o*b_o
	 */
	mq_poly_add_mulX_inner(d, d, t2, ls + 1, logn - 1);
}

/*
 * TODO: make a specialized squaring function, which could be faster
 * than the plain multiplication routine.
 */
static void
mq_poly_sqr_inner(uint8_t *d, const uint8_t *a,
	int ls, unsigned logn, uint8_t *tmp)
{
	mq_poly_mul_inner(d, a, a, ls, logn, tmp);
}

/*
 * Polynomial inversion: we split a into even and odd coefficients:
 *    a = a_e(X^2) + X*a_o(X^2)
 * with a_e and a_o being half-degree. We define:
 *    adj(a) = a_e(X^2) - X*a_o(X^2)
 * Then:
 *    a*adj(a) = a_e^2(X^2) - X^2*a_o^2(X^2)
 *             = (a_e^2 - X*a_o^2)(X^2)
 * which is a half-degree polynomial.
 *
 * Thus:
 *    1/a = adj(a)*(1 / (a*adj(a)))
 * so we reduced the problem of inverting a of degree n into inverting
 * a*adj(a) of degree n/2. We just apply the process recursively
 * until we reach degree 1.
 *
 * Note that 1/(a*adj(a)) is really half-degree, so the multiplication
 * by adj(a) can be done with two half-degree multiplications.
 *
 * Size of tmp[]: 2*n bytes (with n = top-level degree).
 *
 * Return value: 1 on success, 0 on error. On error (a[] is not invertible),
 * contents of d[] are unpredictable.
 */
static int
mq_poly_inv_inner(uint8_t *d, const uint8_t *a, int ls,
	unsigned logn, uint8_t *tmp)
{
	uint8_t *t1, *t2;
	int r;

	if (logn == 1) {
		unsigned a0, a1, x, y;

		a0 = a[0 << ls];
		a1 = a[1 << ls];
		x = a0 * a0 + a1 * a1;

		/*
		 * x is invertible modulo 256 if and only if it is odd.
		 */
		r = (int)(x & 1);

		/*
		 * If x*y = 1 + u*2^k, then:
		 *   x*(y*(2-x*y)) = (1 + u*2^k)*(1 - u*2^k)
		 *                 = 1 - (u^2)*2^(2*k)
		 *                 = 1 mod 2^(2*k)
		 * Inverse of s modulo 4 is itself:
		 *   1*1 = 1 mod 4
		 *   3*3 = 1 mod 4
		 * Thus, we apply the rule above twice on x, to get an
		 * inverse modulo 4*4 = 16, and then modulo 16*16 = 256.
		 */
		y = x;
		y *= 2 - (x * y);
		y *= 2 - (x * y);

		/*
		 * 1/(a0 + X*a1) = (a0 - X*a1) / (a0^2 - X^2*a1^2)
		 *               = (a0 - X*a1) / (a0^2 + a1^2)
		 * (since we work modulo X^2+1)
		 */
		d[0 << ls] = a0 * y;
		d[1 << ls] = -a1 * y;
		return r;
	}

	t1 = tmp + (1 << ls);
	t2 = t1 + ((size_t)1 << (logn + ls));

	/*
	 * t1 <- a*adj(a)
	 */
	mq_poly_sqr_inner(t1, a, ls + 1, logn - 1, tmp);
	mq_poly_sqr_inner(t2, a + (1 << ls), ls + 1, logn - 1, tmp);
	mq_poly_sub_mulX_inner(t1, t1, t2, ls + 1, logn - 1);

	/*
	 * t1 <- 1/a*adj(a)
	 */
	r = mq_poly_inv_inner(t1, t1, ls + 1, logn - 1, tmp);

	/*
	 * d_e <- t1*a_e
	 * d_o <- -t1*a_o
	 */
	mq_poly_mul_inner(d, a, t1, ls + 1, logn - 1, tmp);
	mq_poly_mul_inner(d + (1 << ls),
		a + (1 << ls), t1, ls + 1, logn - 1, tmp);
	mq_poly_neg_inner(d + (1 << ls), d + (1 << ls), ls + 1, logn - 1);

	return r;
}

/*
 * Wrappers for the case of ls = 1 (minimal stride).
 */

static inline void
mq_poly_add(uint8_t *d, const uint8_t *a, const uint8_t *b, unsigned logn)
{
	mq_poly_add_inner(d, a, b, 0, logn);
}

static inline void
mq_poly_sub(uint8_t *d, const uint8_t *a, const uint8_t *b, unsigned logn)
{
	mq_poly_sub_inner(d, a, b, 0, logn);
}

static inline void
mq_poly_mul(uint8_t *d, const uint8_t *a, const uint8_t *b, unsigned logn,
	uint8_t *tmp)
{
	mq_poly_mul_inner(d, a, b, 0, logn, tmp);
}

static inline int
mq_poly_inv(uint8_t *d, const uint8_t *a, unsigned logn, uint8_t *tmp)
{
	return mq_poly_inv_inner(d, a, 0, logn, tmp);
}

/*
 * Multiply polynomial a[] by 1+X+X^2+X^3+...+X^(n-1).
 *
 *   d[0]   = a[0] - a[1] - a[2] - a[3] - ... - a[n - 1]
 *   d[1]   = a[0] + a[1] - a[2] - a[3] - ... - a[n - 1]
 *   d[2]   = a[0] + a[1] + a[2] - a[3] - ... - a[n - 1]
 *     ...
 *   d[n-1] = a[0] + a[1] + a[2] + a[3] + ... + a[n - 1]
 *
 * Thus, d[n - 1] is the sum of all a[i], and:
 *   d[i] = d[i + 1] - 2*a[i + 1]
 * Equivalently:
 *   d[i] = d[i - 1] + 2*a[i]
 * for all i >= 1.
 *
 * This allows efficient computation, in O(n) operations and with no
 * need for extra storage.
 */
static void
mq_poly_mul_ones(uint8_t *d, const uint8_t *a, unsigned logn)
{
	size_t u, n;
	unsigned t;

	n = (size_t)1 << logn;
	t = a[0];
	for (u = 1; u < n; u ++) {
		t -= a[u];
	}
	d[0] = t;
	for (u = 1; u < n; u ++) {
		t += a[u] << 1;
		d[u] = t;
	}
}

/*
 * Multiply a polynomial by a constant.
 */
static void
mq_poly_mulconst(uint8_t *d, const uint8_t *a, unsigned c, unsigned logn)
{
	size_t u, n;

	n = (size_t)1 << logn;
	for (u = 0; u < n; u ++) {
		d[u] = a[u] * c;
	}
}

/* see inner.h */
int
bat_make_public_128(uint8_t *h, const int8_t *f, const int8_t *g,
	unsigned logn, uint32_t *tmp)
{
	size_t u, n;
	uint8_t *t1, *t2;
	int r;

	n = (size_t)1 << logn;
	t1 = (uint8_t *)tmp;
	t2 = t1 + n;

	/*
	 * t1 <- 1/f
	 */
	r = mq_poly_inv(t1, (const uint8_t *)f, logn, t2);

	/*
	 * h <- t1*g = g/f
	 */
	mq_poly_mul(h, t1, (const uint8_t *)g, logn, t2);

	/*
	 * Reduce coefficients modulo 128.
	 */
	for (u = 0; u < n; u ++) {
		h[u] &= 0x7F;
	}
	return r;
}

/* see inner.h */
uint32_t
bat_encrypt_128(int8_t *c, const uint8_t *sbuf,
	const uint8_t *h, unsigned logn, uint32_t *tmp)
{
	size_t u, n;
	uint8_t *t1, *t2;

	n = (size_t)1 << logn;
	t1 = (uint8_t *)tmp;
	t2 = t1 + n;

	/*
	 * Expand sbuf[] into polynomial s (in t1).
	 */
	for (u = 0; u < n; u ++) {
		t1[u] = (sbuf[u >> 3] >> ((unsigned)u & 7)) & 1;
	}

	/*
	 * t1 <- h*s
	 */
	mq_poly_mul(t1, h, t1, logn, t2);

	/*
	 * c = round((h*s) / 2)
	 * Coefficients of h*s must be reduced modulo 128, into -63..+64.
	 * Rounding is toward +inf, thus the result is in -31..+32.
	 */
	for (u = 0; u < n; u ++) {
		c[u] = (((t1[u] + 63) & 0x7F) >> 1) - 31;
	}

	/*
	 * Since coefficients of e' (centered error vector) are in
	 * {-1/2,+1/2}, and those of s' are also in {-1/2,+1/2}, the
	 * norm of vector (gamma*s',e') is always equal to
	 * sqrt((gamma^2 + 1)*(n/4)) = sqrt(n/2), thus always acceptable.
	 */
	return 1;
}

/* see inner.h */
void
bat_decrypt_128(uint8_t *sbuf, const int8_t *c,
	const int8_t *f, const int8_t *g, const int8_t *F, const int8_t *G,
	const int32_t *w, unsigned logn, uint32_t *tmp)
{
	/*
	 * q = 128, Q = 2, q' = 3329, k = 2.
	 *
	 * Decapsulation algorithm:
	 *
	 *   c <- k*c
	 *   c' <- (Q*f*c - f*ones - g*ones) mod q*Q
	 *   c'' <- (q'*Q*F*c - q'*F*ones - q'*G*ones - c'*w) mod q*q'*Q
	 *   e' = (-Gd*c' + g*c'') / (q*q'*Q)
	 *   s' = (Fd*c' - f*c'') / (q*q'*Q)
	 *   e = e' + (1/2)*ones
	 *   s = s' + (1/2)*ones
	 *
	 * We don't need to recompute e, only s. q*Q = 256, which is
	 * natively supported by the code in this file.
	 */
	size_t u, n;
	uint8_t *t1, *t2, *t3, *t4;
	uint16_t *tw1, *tw3, *tw4;

	n = (size_t)1 << logn;
	t1 = (uint8_t *)tmp;
	t2 = t1 + n;
	t3 = t2 + n;
	t4 = t3 + n;

	tw1 = (uint16_t *)t1;
	tw3 = (uint16_t *)t3;
	tw4 = tw3 + n;

	/*
	 * c <- k*c (implicit)
	 * t2 <- Q*c
	 */
	for (u = 0; u < n; u ++) {
		t2[u] = (uint8_t)(c[u] << 2);
	}

	/*
	 * t1 <- c' = (Q*f*c - f*ones - g*ones) mod q*Q
	 */
	mq_poly_mul(t1, t2, (const uint8_t *)f, logn, t3);
	mq_poly_add(t3, (const uint8_t *)f, (const uint8_t *)g, logn);
	mq_poly_mul_ones(t3, t3, logn);
	mq_poly_sub(t1, t1, t3, logn);

	/*
	 * t2 <- (q'*Q*F*c - q'*F*ones - q'*G*ones - c'*w) mod q*Q
	 */
	mq_poly_mul(t2, t2, (const uint8_t *)F, logn, t4);
	mq_poly_add(t3, (const uint8_t *)F, (const uint8_t *)G, logn);
	mq_poly_mul_ones(t3, t3, logn);
	mq_poly_sub(t2, t2, t3, logn);
	mq_poly_mulconst(t2, t2, 64513 & 0xFF, logn);
	for (u = 0; u < n; u ++) {
		t3[u] = (uint8_t)w[u];
	}
	mq_poly_mul(t3, t1, t3, logn, t4);
	mq_poly_sub(t2, t2, t3, logn);

	/*
	 * tw3 <- -c'*w mod q'
	 * This involves rebulding c' mod q' in tw3 first.
	 */
	for (u = 0; u < n; u ++) {
		*(int16_t *)&tw3[u] = (int)((t1[u] + 127) & 0xFF) - 127;
	}
	bat_polyqp_mulneg((int16_t *)tw3, (int16_t *)tw3, w,
		logn, (uint32_t *)tw4);

	/*
	 * At that point, we have:
	 *   t1:  c' mod q*Q   (0..255)
	 *   t2:  c'' mod q*Q  (0..255)
	 *   tw3: c'' mod q'   (-32256..+32256)
	 *
	 * We now want to assemble c'' mod 257 in tw3. We do so by
	 * applying the CRT between t2 and tw3.
	 */
	for (u = 0; u < n; u ++) {
		uint16_t z;
		uint32_t x0, x1, x;
		int32_t y;

		/*
		 * We ensure that we get a positive value by adding
		 * 64513 to the coefficient from c''.
		 */
		x0 = t2[u];
		z = tw3[u];
		x1 = (uint32_t)(*(int16_t *)&z + 64513);

		/*
		 * CRT reconstruction: If:
		 *   x = x0 mod q*Q
		 *   x = x1 mod q'
		 * then:
		 *   x = ((1/q') * (x0 - x1) mod q*Q) * q' + x1
		 * We have q*Q = 256; since 64513 = 252*256 + 1, the value
		 * 1/q' mod q*Q is trivial.
		 */
		x = x1 + ((x0 - x1) & 0xFF) * 64513;

		/*
		 * x is in 0..16515327; we should normalize it to
		 * -8257663..+8257664.
		 */
		y = (int32_t)x - (int32_t)(16515328
			& -((uint32_t)(8257664 - x) >> 31));

		/*
		 * For reduction modulo 257, we ensure a positive value by
		 * adding 32131*257 = 8257667.
		 */
		tw3[u] = m257_tomonty((uint32_t)(y + 8257667));
	}

	/*
	 * We assemble c' mod 257 in tw4 (then moved to tw1).
	 */
	for (u = 0; u < n; u ++) {
		/*
		 * For x in 0..255, normalization to -127..+128 is done
		 * by computing ((x + 127) % 256) - 127. But we then want
		 * to add 257 to get a positive value for reduction modulo
		 * 257.
		 */
		tw4[u] = m257_tomonty(((t1[u] + 127) & 0xFF) + 130);
	}
	memcpy(tw1, tw4, n * sizeof *tw4);

	/*
	 * We have c' mod 257 in tw1, and c'' mod 257 in tw3. We
	 * use the mod 257 code to obtain q*q'*Q*s'.
	 */
	bat_finish_decapsulate_257(tw1, tw3, f, F, w, logn, (uint32_t *)tw4);

	/*
	 * If the ciphertext is correct and the decapsulation worked well,
	 * then s' has coefficients in {-1/2,+1/2} and the coefficients of
	 * s are obtained by adding 1/2. We have the coefficients of
	 * q*q'*Q*s' in tw1[], in Montgomery representation modulo 257:
	 *
	 *    s'   s   tw1[]
	 *  -1/2   0      3     (-q*q'*Q/2 = 3 mod 257)
	 *  +1/2   1    254     (+q*q'*Q/2 = 254 mod 257)
	 *
	 * Thus, we just need to look at the least significant bit of each
	 * value in tw1[] to get the coefficients of s.
	 */
	memset(sbuf, 0, (n + 7) >> 3);
	for (u = 0; u < n; u ++) {
		sbuf[u >> 3] |= (1 - (tw1[u] & 1)) << (u & 7);
	}
}

/* see inner.h */
int
bat_rebuild_G_128(int8_t *G,
	const int8_t *f, const int8_t *g, const int8_t *F,
	unsigned logn, uint32_t *tmp)
{
	size_t u, n;
	uint8_t *t1, *t2, *t3;
	int lim;

	n = (size_t)1 << logn;
	t1 = (uint8_t *)tmp;
	t2 = t1 + n;
	t3 = t2 + n;

	/*
	 * We have g*F - f*G = q; therefore: G = (g*F - q) / f
	 *
	 * We compute modulo 256; note that if f is invertible modulo
	 * 128, it will be invertible modulo 256, and vice versa.
	 */
	if (!mq_poly_inv(t1, (const uint8_t *)f, logn, t3)) {
		return 0;
	}
	mq_poly_mul(t2, (const uint8_t *)g, (const uint8_t *)F, logn, t3);
	t2[0] -= 128;
	mq_poly_mul(t1, t1, t2, logn, t3);

	/*
	 * Normalize coefficients of G around 0, and check that they
	 * are within the expected bounds.
	 */
	lim = (1 << (bat_max_FG_bits[logn] - 1)) - 1;
	for (u = 0; u < n; u ++) {
		int x;

		x = ((t1[u] + 127) & 0xFF) - 127;
		if (x < -lim || x > +lim) {
			return 0;
		}
		G[u] = (int8_t)x;
	}
	return 1;
}

/*
 * Values modulo q are in 0..127, which naturally encodes over exactly
 * 7 bits.
 */

/* see inner.h */
size_t
bat_encode_128(void *out, size_t max_out_len,
	const uint8_t *x, unsigned logn)
{
	size_t u, v, n, out_len;
	uint8_t *buf;

	n = (size_t)1 << logn;
	out_len = ((7 * n) + 7) >> 3;
	if (out == NULL) {
		return out_len;
	}
	if (max_out_len < out_len) {
		return 0;
	}
	buf = out;
	if (n == 2) {
		uint32_t w;

		w = (uint32_t)x[0]
			| ((uint32_t)x[1] << 7);
		enc16le(buf, w);
		return 2;
	} else if (n == 4) {
		uint32_t w;

		w = (uint32_t)x[0]
			| ((uint32_t)x[1] << 7)
			| ((uint32_t)x[2] << 14)
			| ((uint32_t)x[3] << 21);
		enc32le(buf, w);
		return 4;
	} else {
		v = 0;
		for (u = 0; (u + 8) <= n; u += 8) {
			uint32_t w0, w1;

			w0 = (uint32_t)x[u]
				| ((uint32_t)x[u + 1] << 7)
				| ((uint32_t)x[u + 2] << 14)
				| ((uint32_t)x[u + 3] << 21);
			w1 = (uint32_t)x[u + 4]
				| ((uint32_t)x[u + 5] << 7)
				| ((uint32_t)x[u + 6] << 14)
				| ((uint32_t)x[u + 7] << 21);
			enc32le(buf + v, w0 | (w1 << 28));
			enc24le(buf + v + 4, w1 >> 4);
			v += 7;
		}
		return v;
	}
}

/* see inner.h */
size_t
bat_decode_128(uint8_t *x, unsigned logn,
	const void *in, size_t max_in_len)
{
	size_t u, v, n, in_len;
	const uint8_t *buf;
	uint32_t r;

	n = (size_t)1 << logn;
	in_len = ((7 * n) + 7) >> 3;
	if (max_in_len < in_len) {
		return 0;
	}
	buf = in;
	if (n == 2) {
		uint32_t w;

		w = dec16le(buf);
		x[0] = w & 0x7F;
		x[1] = (w >> 7) & 0x7F;
		r = (uint32_t)((w >> 14) - 1) >> 31;
		v = 2;
	} else if (n == 4) {
		uint32_t w;

		w = dec32le(buf);
		x[0] = w & 0x7F;
		x[1] = (w >> 7) & 0x7F;
		x[2] = (w >> 14) & 0x7F;
		x[3] = (w >> 21) & 0x7F;
		r = (uint32_t)((w >> 28) - 1) >> 31;
		v = 4;
	} else {
		v = 0;
		for (u = 0; (u + 8) <= n; u += 8) {
			uint32_t w0, w1;

			w0 = dec32le(buf + v);
			w1 = dec24le(buf + v + 4);
			v += 7;
			x[u + 0] = w0 & 0x7F;
			x[u + 1] = (w0 >> 7) & 0x7F;
			x[u + 2] = (w0 >> 14) & 0x7F;
			x[u + 3] = (w0 >> 21) & 0x7F;
			x[u + 4] = ((w0 >> 28) | (w1 << 4)) & 0x7F;
			x[u + 5] = (w1 >> 3) & 0x7F;
			x[u + 6] = (w1 >> 10) & 0x7F;
			x[u + 7] = (w1 >> 17) & 0x7F;
		}
		r = 1;
	}
	return v & -(size_t)r;
}

/*
 * Ciphertext values are in -31..+32 range; for value v, we encode v+31
 * over 6 bits (v+31 is in 0..63).
 */

/* see inner.h */
size_t
bat_encode_ciphertext_128(void *out, size_t max_out_len,
	const int8_t *c, unsigned logn)
{
	size_t u, v, n, out_len;
	uint8_t *buf;

	n = (size_t)1 << logn;
	out_len = ((6 * n) + 7) >> 3;
	if (out == NULL) {
		return out_len;
	}
	if (max_out_len < out_len) {
		return 0;
	}
	buf = out;
	v = 0;
	for (u = 0; (u + 4) <= n; u += 4) {
		uint32_t w;

		w = (uint32_t)(c[u] + 31)
			| ((uint32_t)(c[u + 1] + 31) << 6)
			| ((uint32_t)(c[u + 2] + 31) << 12)
			| ((uint32_t)(c[u + 3] + 31) << 18);
		enc24le(buf + v, w);
		v += 3;
	}
	if (u < n) {
		uint32_t w;

		w = (uint32_t)(c[u] + 31)
			| ((uint32_t)(c[u + 1] + 31) << 6);
		enc16le(buf + v, w);
		v += 2;
	}
	return v;
}

/* see inner.h */
size_t
bat_decode_ciphertext_128(int8_t *c, unsigned logn,
	const void *in, size_t max_in_len)
{
	size_t u, v, n, in_len;
	const uint8_t *buf;
	uint32_t r;

	n = (size_t)1 << logn;
	in_len = ((6 * n) + 7) >> 3;
	if (max_in_len < in_len) {
		return 0;
	}
	buf = in;
	v = 0;
	r = 1;
	for (u = 0; (u + 4) <= n; u += 4) {
		uint32_t w;

		w = dec24le(buf + v);
		v += 3;
		c[u + 0] = (int)(w & 0x3F) - 31;
		c[u + 1] = (int)((w >> 6) & 0x3F) - 31;
		c[u + 2] = (int)((w >> 12) & 0x3F) - 31;
		c[u + 3] = (int)((w >> 18) & 0x3F) - 31;
	}
	if (u < n) {
		uint32_t w;

		w = dec16le(buf + v);
		v += 2;
		c[u + 0] = (int)(w & 0x3F) - 31;
		c[u + 1] = (int)((w >> 6) & 0x3F) - 31;
		r &= (uint32_t)((w >> 12) - 1) >> 31;
	}
	return v & -(size_t)r;
}
