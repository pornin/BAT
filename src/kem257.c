#include "inner.h"

/*
 * We use the mod 257 code here.
 */
#define Q   257
#include "modgen.c"

/* see inner.h */
int
bat_make_public_257(uint16_t *h, const int8_t *f, const int8_t *g,
	unsigned logn, uint32_t *tmp)
{
	size_t u, n;
	uint16_t *ft, *gt;

	n = (size_t)1 << logn;
	ft = (uint16_t *)tmp;
	gt = ft + n;
	for (u = 0; u < n; u ++) {
		ft[u] = mq_set(f[u]);
		gt[u] = mq_set(g[u]);
	}
	NTT(ft, ft, logn);
	if (!mq_poly_inv_ntt(ft, ft, logn)) {
		return 0;
	}
	NTT(gt, gt, logn);
	mq_poly_mul_ntt(ft, ft, gt, logn);
	iNTT(ft, ft, logn);
	for (u = 0; u < n; u ++) {
		h[u] = mq_unorm(ft[u]);
	}
	return 1;
}

/* see inner.h */
uint32_t
bat_encrypt_257(int8_t *c, const uint8_t *sbuf,
	const uint16_t *h, unsigned logn, uint32_t *tmp)
{
	size_t u, n;
	uint16_t *t1, *t2;

	n = (size_t)1 << logn;
	t1 = (uint16_t *)tmp;
	t2 = t1 + n;

	/*
	 * Get h into t2, in NTT representation.
	 */
	for (u = 0; u < n; u ++) {
		t2[u] = mq_set(h[u]);
	}
	NTT(t2, t2, logn);

	/*
	 * Coefficients of polynomial s are {0,1}, extracted from sbuf[].
	 */
	for (u = 0; u < n; u ++) {
		t1[u] = mq_set((sbuf[u >> 3] >> (u & 7)) & 1);
	}

	/*
	 * t1 <- h*s
	 * We have the NTT representation of h in t2.
	 */
	NTT(t1, t1, logn);
	mq_poly_mul_ntt(t1, t1, t2, logn);
	iNTT(t1, t1, logn);

	/*
	 * c <- round((h*s mod q)/k)
	 * (for q = 257, we have k = 2).
	 *
	 * Rounding is toward +infty, so that error e = k*c - (h*s mod q)
	 * has coefficients in {0, 1} only.
	 */
	for (u = 0; u < n; u ++) {
		c[u] = ((mq_snorm(t1[u]) + 129) >> 1) - 64;
	}

	/*
	 * Ciphertext is acceptable if and only if the norm
	 * of (gamma*s',e') is not greater than
	 * 1.08*sqrt((n/2)*gamma^2). Note that e' and s' are
	 * centered on 0, i.e. e' = e - E(e) and s' = s - E(s).
	 * With k = 2, we have gamma = 1, and all coefficients
	 * of e' and s' can only be 1/2 or -1/2. Therefore,
	 * the norm of (gamma*s',e') can only be sqrt(n/2),
	 * therefore always below the bound.
	 */
	return 1;
}

/* see inner.h */
void
bat_decrypt_257(uint8_t *sbuf, const int8_t *c,
	const int8_t *f, const int8_t *g, const int8_t *F, const int8_t *G,
	const int32_t *w, unsigned logn, uint32_t *tmp)
{
	/*
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
	 * If the ciphertext is correct, then it must be that all
	 * coefficients of s are in {0,1}, and all coefficients of e are
	 * in -((k/2)-1)..(k/2) (for k even).
	 *
	 * When q = 257, we have Q = 2 and k = 2. This simplifies some
	 * things:
	 *
	 *   - f*ones mod 2  is a constant polynomial: all its
	 *     coefficients are equal to the "parity" of f, i.e. the sum
	 *     of its coefficients mod 2. The keygen implementation
	 *     produces only polynomials f and g with odd parity, but an
	 *     externally provided key may have f or g with an even
	 *     parity (not both can be even, though, otherwise the NTRU
	 *     equation g*F - f*G = q has no solution).
	 *
	 *   - We compute c' modulo q, and modulo Q = 2. Since f*ones
	 *     and g*ones are constant polynomials modulo 2, and
	 *     Q*f*c = 0 mod Q, c' mod Q is a constant polynomial (all
	 *     coefficients are 0, or all coefficients are 1). We
	 *     get c' mod q, and adjust each coefficient in order to have
	 *     the correct parity.
	 *
	 *   - c'' is computed modulo q, q' and Q. Modulo Q = 2:
	 *        q'*Q*F*c = 0 mod Q
	 *        q'*F*ones mod Q  is a constant polynomial
	 *        q'*G*ones mod Q  is a constant polynomial
	 *        c' mod Q  is a constant polynomial, either 0 or ones;
	 *        thus, c'*w is a constant polynomial modulo Q.
	 *
	 *     Modulo q', q'*Q*F*c, q'*F*ones and q'*G*ones are zero,
	 *     so only c'*w has to be computed modulo q'.
	 *
	 *     Modulo q, all elements must be computed.
	 *
	 *     When we have c'' modulo q, q' and Q, we use the CRT to get
	 *     the value modulo q*q'*Q.
	 *
	 * Once we have c' and c'', we can compute 2*q*q'*Q*e' and
	 * 2*q*q'*Q*s'. These are nominally plain integers, but we do
	 * not need the full values; we just want to distinguish the
	 * possible values for the coefficients of s' and e' (these
	 * coefficients can only be 0 or 1). Thus, we can do the
	 * computation modulo any prime p which is not 2, q or q'; in
	 * practice, we use the "other q" (i.e. 769) since we already
	 * have the code for computations modulo that prime.
	 */
	size_t u, n;
	uint16_t *t1, *t2, *t3, *t4;
	unsigned par_fg, par_FG, par_w, cp2, cs2;

	n = (size_t)1 << logn;

	t1 = (uint16_t *)tmp;
	t2 = t1 + n;
	t3 = t2 + n;
	t4 = t3 + n;

	/*
	 * In the BAT specification, algorithm 3.2 ("Decode") expects
	 * as polynomial 'c' what is really 'k*c' in algorithm 4.3
	 * ("Decapsulate"). With q = 257, we have k = 2. Moreover, we
	 * want Q*c, with Q = 2; thus, we compute 4*c here.
	 */
	for (u = 0; u < n; u ++) {
		t1[u] = mq_set(4 * c[u]);
	}

	/*
	 * We have Q*c in t1; convert it to NTT (mod q).
	 */
	NTT(t1, t1, logn);

	/*
	 * Get NTT representations of f and f+g, in t2 and t3, respectively.
	 */
	for (u = 0; u < n; u ++) {
		t2[u] = mq_set(f[u]);
		t3[u] = mq_set(f[u] + g[u]);
	}
	NTT(t2, t2, logn);
	NTT(t3, t3, logn);

	/*
	 * t2 <- Q*f*c mod q  (NTT)
	 */
	mq_poly_mul_ntt(t2, t1, t2, logn);

	/*
	 * t2 <- c' mod q  (NTT)
	 * t3 is scratch
	 */
	mq_poly_mul_ones_ntt(t3, t3, logn);
	mq_poly_sub(t2, t2, t3, logn);

	/*
	 * t1 <- q'*Q*F*c mod q  (NTT)
	 * t3 <- q'*F mod q      (NTT)
	 */
	for (u = 0; u < n; u ++) {
		t3[u] = mq_set((int32_t)F[u] * (64513 % 257));
	}
	NTT(t3, t3, logn);
	mq_poly_mul_ntt(t1, t1, t3, logn);

	/*
	 * t1 <- q'*Q*F*c - q'*F*ones mod q  (NTT)
	 * t3 is scratch
	 */
	mq_poly_mul_ones_ntt(t3, t3, logn);
	mq_poly_sub(t1, t1, t3, logn);

	/*
	 * t1 <- q'*Q*F*c - q'*F*ones - q'*G*ones mod q  (NTT)
	 */
	for (u = 0; u < n; u ++) {
		t3[u] = mq_set((int32_t)G[u] * (64513 % 257));
	}
	NTT(t3, t3, logn);
	mq_poly_mul_ones_ntt(t3, t3, logn);
	mq_poly_sub(t1, t1, t3, logn);

	/*
	 * t1 <- q'*Q*F*c - q'*F*ones - q'*G*ones - c'*w mod q  (NTT)
	 */
	for (u = 0; u < n; u ++) {
		t3[u] = mq_set(w[u]);
	}
	NTT(t3, t3, logn);
	mq_poly_mul_ntt(t3, t2, t3, logn);
	mq_poly_sub(t1, t1, t3, logn);

	/*
	 * We won't need NTT representations for c' and c'' any more, so
	 * we convert back to normal (but still Montgomery).
	 */
	iNTT(t1, t1, logn);
	iNTT(t2, t2, logn);

	/*
	 * We now have c'' and c' mod q, in t1 and t2, respectively (both
	 * are in normal representation).
	 *
	 * The actual values (with plain integer coefficients) of c'
	 * and c'' are the results of modular reduction and normalization
	 * of c' and c'' modulo q*Q and q*q'*Q, respectively. We must
	 * thus get c' modulo Q, and c'' modulo Q and modulo q', in order
	 * to access the actual coefficients.
	 *
	 *
	 * Since Q = 2, there are shortcuts for computing modulo Q. The
	 * product of any polynomial by the 'ones' polynomial yields
	 * either zero, or the 'ones' polynomial, depending on the
	 * "parity" of the source polynomial (the sum of its
	 * coefficients modulo 2). Therefore, modulo Q:
	 *
	 *   c' = Q*f*c - (f+g)*ones mod 2 = parity(f+g)*ones modulo 2
	 *   c'' = q'*Q*F*c - q'*(F+G)*ones - c'*w mod 2
	 *       = (parity(F+G) + parity(f+g)*parity(w) mod 2)*ones
	 *
	 * Therefore, we only need to get the parities of f, g, F, G and
	 * w to get c' and c'' modulo Q.
	 *
	 * Note that our keygen algorithm only produces f and g such
	 * that parity(f) = 1 and parity(g) = 1. Thus, for our private
	 * keys, we have parity(f+g) = 0. However, this may not be the
	 * case for other, externally generated private keys (the NTRU
	 * equation g*F - f*G = q cannot be solved for an odd q if f and
	 * g both have parity 0, but it is possible to have one of f or
	 * g to be of parity 0). Therefore, we apply the complete
	 * treatment below (the extra cost is negligible in practice).
	 */
	par_fg = 0;
	par_FG = 0;
	par_w = 0;
	for (u = 0; u < n; u ++) {
		par_fg += (unsigned)f[u] + (unsigned)g[u];
		par_FG += (unsigned)F[u] + (unsigned)G[u];
		par_w += (unsigned)w[u];
	}
	par_fg &= 1;
	par_FG &= 1;
	par_w &= 1;

	cp2 = par_fg;
	cs2 = par_FG ^ (par_fg & par_w);

	/*
	 * Set t2 to coefficients of c'. We have c' mod q in t2, so we
	 * just need to adjust them to account for the parity we just
	 * computed in cp2:
	 *  - coefficient must be normalized to -128..+128
	 *  - if the parity is wrong, we must add 257 (if the value is
	 *    negative or zero) or subtract 257 (if the value is strictly
	 *    positive), so that the result is in -256..+257.
	 */
	for (u = 0; u < n; u ++) {
		uint32_t x;

		/*
		 * Get next coefficient of c' into x, normalized in
		 * -128..+128.
		 */
		x = (uint32_t)mq_snorm(t2[u]);

		/*
		 * Adjust x to get the right value modulo 2: if it has
		 * the wrong parity, then we must add or subtract 257.
		 * We add 257 if the value is negative or zero, subtract
		 * 257 if it is negative.
		 */
		x += -(uint32_t)((x ^ cp2) & 1u)
			& -(uint32_t)257
			& (((x - 1) >> 16) & (2 * 257));

		/*
		 * We write back x as an unsigned 16-bit integer, but with
		 * two's complement notation for negative values, which will
		 * be fine since uint16_t and int16_t have compatible
		 * memory layouts, and bat_polyqp_mulneg() expects int16_t
		 * values.
		 */
		t2[u] = (uint16_t)x;
	}

	/*
	 * We want:
	 *   c'' = -c'*w mod q'
	 * We compute that value into t3.
	 */
	bat_polyqp_mulneg((int16_t *)t3, (int16_t *)t2, w, logn,
		(uint32_t *)t4);

	/*
	 * At that point:
	 *    t1    c'' mod q  (normal representation, Montgomery)
	 *    t2    c'         (normal representation, signed)
	 *    t3    c'' mod q' (normal representation, signed)
	 *    cs2   c'' mod 2  (to multiply by 'ones')
	 * We now compute c'' itself in normal representation and signed
	 * integers by combining the coefficients modulo q, q' and 2,
	 * stored in t1[], t3[] and cs2, respectively. This uses the CRT.
	 */
	for (u = 0; u < n; u ++) {
		uint32_t y0, y1, x;

		/*
		 * If:
		 *    y = y0 mod q
		 *    y = y1 mod q'
		 * Then:
		 *    y = ((1/q') * (y0 - y1) mod q) * q' + y1
		 * We need y0 and y1 in 0..q-1 and 0..q'-1, respectively.
		 */
		y0 = mq_unorm(t1[u]);
		y1 = (uint32_t)*(int16_t *)&t3[u];
		y1 += 64513 & (y1 >> 16);

		/*
		 * The Montgomery representation of 1/q' mod q is 43
		 * (with q = 257 and q' = 64513). We add 64764 = 252*257 to
		 * ensure that the value provided to mq_montyred() is in
		 * the proper range (max value will be 43*(256+64764)).
		 */
		x = mq_montyred(43 * (64764 + y0 - y1));

		/*
		 * Value x is in 1..q range. We need to normalize value
		 * q to 0.
		 */
		x &= (uint32_t)(x - Q) >> 16;

		/*
		 * Compute value modulo q*q', in 0..q*q'-1 range.
		 */
		x = (x * 64513) + (uint32_t)y1;

		/*
		 * If x = 0, set it to q*q' = 16579841.
		 */
		x += 16579841 & -((uint32_t)(x - 1) >> 31);

		/*
		 * Adjust parity to get value modulo 2*q*q': we subtract
		 * q*q' = 16579841 if the value has the wrong parity.
		 */
		x -= 16579841 & -(uint32_t)((x & 1) ^ cs2);

		/*
		 * Value is now in -8289920..+8289920, which is the correct
		 * normalized range. Since we will hand it over to the
		 * module that computes modulo 769, we pre-reduce it
		 * modulo 769. We ensure a positive value by adding
		 * 8290589 = 769 * 10781.
		 */
		t1[u] = m769_tomonty(x + 8290589);
	}

	/*
	 * We now have c' and c'' in t2 and t1, respectively. c' uses
	 * signed integers; c'' is in Montgomery representation modulo 769.
	 * We convert c' the Montgomery representation modulo 769 as
	 * well.
	 */
	for (u = 0; u < n; u ++) {
		t2[u] = m769_tomonty(*(int16_t *)&t2[u] + 769);
	}

	/*
	 * We need to compute:
	 *   Fd = q'*F - f*w
	 *   s' = 1/(q*q'*Q) (Fd*c' - f*c'')
	 *   s = s' + (1/2)*ones
	 *
	 * We use the mod 769 code to obtain q*q'*Q*s'.
	 */
	bat_finish_decapsulate_769(t2, t1, f, F, w, logn, (uint32_t *)t3);

	/*
	 * If the ciphertext is correct and the decapsulation worked well,
	 * then s' has coefficients in {-1/2,1/2} and the coefficients of
	 * s are obtained by adding 1/2. We have the coefficients of
	 * q*q'*Q*s' in t2[], in Montgomery representation modulo 769:
	 *
	 *    s'   s   t2[]
	 *  -1/2   0    26    (-q*q'*Q/2 = 568 mod 769)
	 *  +1/2   1   743    (+q*q'*Q/2 = 201 mod 769)
	 *
	 * Therefore, we just need to look at the least significant bit
	 * of each value in t2[] to get the coefficients of s.
	 */
	memset(sbuf, 0, (n + 7) >> 3);
	for (u = 0; u < n; u ++) {
		sbuf[u >> 3] |= (t2[u] & 1) << (u & 7);
	}
}

/* see inner.h */
void
bat_finish_decapsulate_257(uint16_t *cp, uint16_t *cs,
	const int8_t *f, const int8_t *F, const int32_t *w, unsigned logn,
	uint32_t *tmp)
{
	/*
	 * Formulas:
	 *   Fd = q'*F - f*w
	 *   q*q'*Q*s = Fd*c' - f*c''
	 * We q' = 3329.
	 */
	size_t u, n;
	uint16_t *t1, *t2;

	n = (size_t)1 << logn;
	t1 = (uint16_t *)tmp;
	t2 = t1 + n;

	/*
	 * Convert c' and c'' to NTT.
	 */
	NTT(cp, cp, logn);
	NTT(cs, cs, logn);

	/*
	 * Load f in t1 and w in t2 and convert both to NTT.
	 */
	for (u = 0; u < n; u ++) {
		t1[u] = mq_set(f[u]);
		t2[u] = mq_set(w[u]);
	}
	NTT(t1, t1, logn);
	NTT(t2, t2, logn);

	/*
	 * cs <- f*c''  (NTT)
	 * t1 <- f*w    (NTT)
	 */
	mq_poly_mul_ntt(cs, cs, t1, logn);
	mq_poly_mul_ntt(t1, t1, t2, logn);

	/*
	 * t2 <- q'*F  (NTT)
	 */
	for (u = 0; u < n; u ++) {
		/*
		 * 771 = 257 * 3.
		 * This addition ensures that (q' mod q)*F[u] becomes a
		 * positive integer, thus in range for mq_tomonty().
		 */
		t2[u] = mq_tomonty((int32_t)F[u] * (64513 % 257) + 771);
	}
	NTT(t2, t2, logn);

	/*
	 * t1 <- Fd = q'*F - f*w  (NTT)
	 */
	mq_poly_sub(t1, t2, t1, logn);

	/*
	 * cp <- Fd*c' - f*c''    (normal representation)
	 */
	mq_poly_mul_ntt(t1, t1, cp, logn);
	mq_poly_sub(t1, t1, cs, logn);
	iNTT(cp, t1, logn);
}

/* see inner.h */
int
bat_rebuild_G_257(int8_t *G,
	const int8_t *f, const int8_t *g, const int8_t *F,
	unsigned logn, uint32_t *tmp)
{
	size_t u, n;
	uint16_t *t1, *t2;
	int lim;

	n = (size_t)1 << logn;
	t1 = (uint16_t *)tmp;
	t2 = t1 + n;

	/*
	 * Load g in t1 and F in t2.
	 */
	for (u = 0; u < n; u ++) {
		t1[u] = mq_set(g[u]);
		t2[u] = mq_set(F[u]);
	}

	/*
	 * Get g*F into t1, in NTT representation.
	 */
	NTT(t1, t1, logn);
	NTT(t2, t2, logn);
	mq_poly_mul_ntt(t1, t1, t2, logn);

	/*
	 * Load f in t2.
	 */
	for (u = 0; u < n; u ++) {
		t2[u] = mq_set(f[u]);
	}

	/*
	 * Compute G in t1 (NTT).
	 */
	NTT(t2, t2, logn);
	if (!mq_poly_inv_ntt(t2, t2, logn)) {
		return 0;
	}
	mq_poly_mul_ntt(t1, t1, t2, logn);

	/*
	 * Get back the coefficients of G. The coefficients are verified
	 * to remain within the expected bound.
	 */
	iNTT(t1, t1, logn);
	lim = (1 << (bat_max_FG_bits[logn] - 1)) - 1;
	for (u = 0; u < n; u ++) {
		int x;

		x = mq_snorm(t1[u]);
		if (x < -lim || x > +lim) {
			return 0;
		}
		G[u] = (int8_t)x;
	}
	return 1;
}

/*
 * Encoding format: 8 values modulo 257 are encoded over 65 bits.
 * Value x_i is split into its 4 low-order bits (xl_i, value in 0..15)
 * and its 5 high-order bits (xh_i, value in 0..16). First 32 bits
 * are the xl_i; then remaining 33 bits are the xh_i (digits in base-17).
 *
 * Eight 65-bit blocks are encoded by grouping the high bits into a common
 * final byte.
 *
 * For small degrees:
 *
 *  - If degree is 2, then encoding is over three bytes, 10 bits per value.
 *  - If degree is 4, then encoding is over five bytes, 10 bits per value.
 *  - If degree is 8, 16 or 32, then one, two or four 65-bit blocks are
 *    used; the ignored bits in the last byte are set to zero.
 */

/*
 * Encode eight values x[0]..x[7] into 65 bits (bits 0..31 are returned
 * in *lo, bits 32..63 in *hi, bit 64 is the function return value).
 */
static inline uint32_t
encode8(uint32_t *lo, uint32_t *hi, const uint16_t *x)
{
	int i;
	uint32_t wl, wh, r, t;

	wl = 0;
	for (i = 0; i < 8; i ++) {
		wl |= (uint32_t)(x[i] & 0x0F) << (4 * i);
	}
	*lo = wl;
	wh = 0;
	for (i = 7; i > 0; i --) {
		wh = (wh * 17) + (uint32_t)(x[i] >> 4);
	}

	/*
	 * Final iteration is:
	 *   *hi = wh*17 + (x[0] >> 4)
	 * which can be written as:
	 *   *hi = wh*16 + wh + (x[0] >> 4)
	 * There will be a carry into the top bit (bit 64 of the final
	 * output) if and only if one of the following holds:
	 *  - wh >= 2^28
	 *  - adding wh + (x[0] >> 4) triggers a carry
	 * Note that wh + (x[0] >> 4) <= 410338688 < 2^29.
	 */

	r = wh >> 28;
	t = wh + (x[0] >> 4);
	wh = (wh << 4) + t;
	r |= ((uint32_t)(wh - t) & ~wh) >> 31;
	*hi = wh;
	return r;
}

#if BAT_AVX2

/*
 * Eight encode8() in parallel. 64 input values (x[0..63]) are encoded
 * into 64 consecutive bytes (buf[0..63]) and the eight extra bits are
 * returned as a packed byte value.
 */
TARGET_AVX2
static inline uint32_t
encode8x8(uint8_t *buf, const uint16_t *x)
{
	__m256i xm4, xm16, xmlo32, xs, x83521;
	__m256i xi0, xi1, xi2, xi3, xa0, xa1, xa2, xa3, xt0, xt1, xt2, xt3;

	xm4 = _mm256_set1_epi16(0x0F);
	xm16 = _mm256_set1_epi32(0xFFFF);
	xmlo32 = _mm256_setr_epi32(-1, 0, 0, 0, -1, 0, 0, 0);
	x83521 = _mm256_set1_epi64x(83521);

	/*
	 * Load the 64 values to encode.
	 */
	xi0 = _mm256_loadu_si256((const void *)(x +  0));
	xi1 = _mm256_loadu_si256((const void *)(x + 16));
	xi2 = _mm256_loadu_si256((const void *)(x + 32));
	xi3 = _mm256_loadu_si256((const void *)(x + 48));

	/*
	 * Repack low nibbles to get low halves.
	 */

	/* Extract low 4 bits of each value. */
	xa0 = _mm256_and_si256(xi0, xm4);
	xa1 = _mm256_and_si256(xi1, xm4);
	xa2 = _mm256_and_si256(xi2, xm4);
	xa3 = _mm256_and_si256(xi3, xm4);

	/* Pack two nibbles into one byte, and shuffle bytes to get
	   one 32-bit word per AVX2 lane (two per register). */
	xs = _mm256_setr_epi8(
		0, 4, 8, 12, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1,
		0, 4, 8, 12, 1, 1, 1, 1,
		1, 1, 1, 1, 1, 1, 1, 1);
	xa0 = _mm256_shuffle_epi8(
		_mm256_or_si256(xa0, _mm256_srli_epi32(xa0, 12)), xs);
	xa1 = _mm256_shuffle_epi8(
		_mm256_or_si256(xa1, _mm256_srli_epi32(xa1, 12)), xs);
	xa2 = _mm256_shuffle_epi8(
		_mm256_or_si256(xa2, _mm256_srli_epi32(xa2, 12)), xs);
	xa3 = _mm256_shuffle_epi8(
		_mm256_or_si256(xa3, _mm256_srli_epi32(xa3, 12)), xs);

	/*
	 * Extract high parts from each input value.
	 */
	xi0 = _mm256_srli_epi16(xi0, 4);
	xi1 = _mm256_srli_epi16(xi1, 4);
	xi2 = _mm256_srli_epi16(xi2, 4);
	xi3 = _mm256_srli_epi16(xi3, 4);

	/*
	 * Combine high parts by pairs.
	 */
	xi0 = _mm256_add_epi32(
		_mm256_add_epi32(xi0, _mm256_srli_epi32(xi0, 16)),
		_mm256_srli_epi32(xi0, 12));
	xi1 = _mm256_add_epi32(
		_mm256_add_epi32(xi1, _mm256_srli_epi32(xi1, 16)),
		_mm256_srli_epi32(xi1, 12));
	xi2 = _mm256_add_epi32(
		_mm256_add_epi32(xi2, _mm256_srli_epi32(xi2, 16)),
		_mm256_srli_epi32(xi2, 12));
	xi3 = _mm256_add_epi32(
		_mm256_add_epi32(xi3, _mm256_srli_epi32(xi3, 16)),
		_mm256_srli_epi32(xi3, 12));

	/*
	 * Clear unusued bits so that we have 32 bits per pair, then
	 * combine pairs. An input pair is x0+17*x1; we need to compute
	 * (x0+17*x1) + 289*(x2+17*x3). The multiplication result may
	 * exceed 16 bits, so we would need to use 32-bit multiplies,
	 * but they have high latencies; instead, we do two "manual"
	 * multiplications by 17.
	 */
	xi0 = _mm256_and_si256(xi0, xm16);
	xi1 = _mm256_and_si256(xi1, xm16);
	xi2 = _mm256_and_si256(xi2, xm16);
	xi3 = _mm256_and_si256(xi3, xm16);

	xt0 = _mm256_add_epi32(xi0, _mm256_slli_epi32(xi0, 4));
	xt1 = _mm256_add_epi32(xi1, _mm256_slli_epi32(xi1, 4));
	xt2 = _mm256_add_epi32(xi2, _mm256_slli_epi32(xi2, 4));
	xt3 = _mm256_add_epi32(xi3, _mm256_slli_epi32(xi3, 4));
	xi0 = _mm256_add_epi32(
		_mm256_add_epi32(xi0, _mm256_srli_epi64(xt0, 32)),
		_mm256_srli_epi64(xt0, 28));
	xi1 = _mm256_add_epi32(
		_mm256_add_epi32(xi1, _mm256_srli_epi64(xt1, 32)),
		_mm256_srli_epi64(xt1, 28));
	xi2 = _mm256_add_epi32(
		_mm256_add_epi32(xi2, _mm256_srli_epi64(xt2, 32)),
		_mm256_srli_epi64(xt2, 28));
	xi3 = _mm256_add_epi32(
		_mm256_add_epi32(xi3, _mm256_srli_epi64(xt3, 32)),
		_mm256_srli_epi64(xt3, 28));

	/*
	 * We have quartets of high parts, four per registers (two
	 * per lane). To combine two quartets, we need to multiply
	 * the high-quartet by 17**4 = 83521. This time, the
	 * mul_epu32 intrinsic is convenient.
	 */
	xi0 = _mm256_add_epi64(
		_mm256_and_si256(xi0, xmlo32),
		_mm256_bsrli_epi128(_mm256_mul_epu32(xi0, x83521), 8));
	xi1 = _mm256_add_epi64(
		_mm256_and_si256(xi1, xmlo32),
		_mm256_bsrli_epi128(_mm256_mul_epu32(xi1, x83521), 8));
	xi2 = _mm256_add_epi64(
		_mm256_and_si256(xi2, xmlo32),
		_mm256_bsrli_epi128(_mm256_mul_epu32(xi2, x83521), 8));
	xi3 = _mm256_add_epi64(
		_mm256_and_si256(xi3, xmlo32),
		_mm256_bsrli_epi128(_mm256_mul_epu32(xi3, x83521), 8));

	/*
	 * At that point, the xa* contain the low 32-bit halves of the
	 * results, and the xi* contain the high 33-bit halves:
	 *
	 *   xa0: lo0 0 0 0 lo1 0 0 0
	 *   xa1: lo2 0 0 0 lo3 0 0 0
	 *   xa2: lo4 0 0 0 lo5 0 0 0
	 *   xa3: lo6 0 0 0 lo7 0 0 0
	 *
	 *   xi0: hi0 t0 0 0 hi1 t1 0 0
	 *   xi1: hi2 t2 0 0 hi3 t3 0 0
	 *   xi2: hi4 t4 0 0 hi5 t5 0 0
	 *   xi3: hi6 t6 0 0 hi7 t7 0 0
	 *
	 * We need to repack everything into:
	 *
	 *   lo0 hi0 lo1 hi1 lo2 hi2 lo3 hi3
	 *   lo4 hi4 lo5 hi5 lo6 hi6 lo7 hi7
	 *   and the t* bits into a byte.
	 */

	/* Pack lo* and hi* together (without the t* bits). */
	xa0 = _mm256_or_si256(xa0, _mm256_slli_epi64(xi0, 32));
	xa1 = _mm256_or_si256(xa1, _mm256_slli_epi64(xi1, 32));
	xa2 = _mm256_or_si256(xa2, _mm256_slli_epi64(xi2, 32));
	xa3 = _mm256_or_si256(xa3, _mm256_slli_epi64(xi3, 32));

	/* Pack into two registers. */
	xa0 = _mm256_or_si256(xa0, _mm256_bslli_epi128(xa1, 8));
	xa2 = _mm256_or_si256(xa2, _mm256_bslli_epi128(xa3, 8));

	/*
	 *   xa0: lo0 hi0 lo2 hi2 lo1 hi1 lo3 hi3
	 *   xa2: lo4 hi4 lo6 hi6 lo5 hi5 lo7 hi7
	 */

	/* Cross-lane permutations to get the right output order. */
	xa0 = _mm256_permute4x64_epi64(xa0, 0xD8);
	xa2 = _mm256_permute4x64_epi64(xa2, 0xD8);

	/*
	 * Write the 64 output bytes.
	 */
	_mm256_storeu_si256((void *)(buf +  0), xa0);
	_mm256_storeu_si256((void *)(buf + 32), xa2);

	/*
	 * Compute the extra byte.
	 */
	xi0 = _mm256_srli_epi64(xi0, 32);
	xi1 = _mm256_srli_epi64(xi1, 32);
	xi2 = _mm256_srli_epi64(xi2, 32);
	xi3 = _mm256_srli_epi64(xi3, 32);
	xi0 = _mm256_or_si256(
		_mm256_or_si256(
			xi0,
			_mm256_slli_epi32(xi1, 2)),
		_mm256_or_si256(
			_mm256_slli_epi32(xi2, 4),
			_mm256_slli_epi32(xi3, 6)));
	xi0 = _mm256_permute4x64_epi64(xi0, 0xD8);
	xi0 = _mm256_or_si256(xi0, _mm256_bsrli_epi128(xi0, 1));
	xi0 = _mm256_or_si256(xi0, _mm256_srli_epi64(xi0, 55));
	return _mm256_cvtsi256_si32(xi0);
}

#endif

/*
 * Decode a 65-bit word (lo is bits 0..31, hi is bits 32..63, tt is bit 64)
 * into eight values. Returned value is 1 on success, 0 on error.
 */
static uint32_t
decode8(uint16_t *x, uint32_t lo, uint32_t hi, uint32_t tt)
{
	/*
	 * Useful facts:
	 *    256 = 1 mod 17
	 *    65536 = 1 mod 17
	 *    2^32 = 1 mod 17
	 *    for any x in 0..69631, (61681 * x) >> 20 = floor(x / 17)
	 *    17 * 4042322161 = 1 mod 2^32
	 */
	uint32_t a;
	int i;

	/*
	 * First iteration: we need to compute the quotient and
	 * remainder of (hi + (tt << 32)) divided by 17.
	 *
	 * We first compute the remainder.
	 */
	a = (hi & 0xFFFF) + (hi >> 16) + tt;   /* 0 <= a <= 131071 */
	a = (a & 0xFF) + (a >> 8);             /* 0 <= a <= 766 */
	a -= 17 * ((61681 * a) >> 20);         /* 0 <= a <= 16 */
	x[0] = (lo & 0x0F) + (a << 4);

	/*
	 * Subtract the remainder, then do the exact division by multiplying
	 * by the inverse of 17 modulo 2^32. We do all computations modulo
	 * 2^32, since we know that the result will be lower than 2^32 (we
	 * don't have to care about the value of tt here, it is implictly
	 * taken into account in the remainder).
	 */
	hi = (hi - a) * 4042322161u;

	/*
	 * Next four iterations are similar, except that we work on a
	 * 32-bit value now (no tt).
	 */
	for (i = 1; i <= 4; i ++) {
		a = (hi & 0xFFFF) + (hi >> 16);   /* 0 <= a <= 131070 */
		a = (a & 0xFF) + (a >> 8);        /* 0 <= a <= 766 */
		a -= 17 * ((61681 * a) >> 20);    /* 0 <= a <= 16 */
		x[i] = ((lo >> (4 * i)) & 0x0F) + (a << 4);
		hi = (hi - a) * 4042322161u;
	}

	/*
	 * At that point, max value for hi is 6049. The next two iterations
	 * can be done with simpler computations.
	 */
	for (i = 5; i <= 6; i ++) {
		a = (61681 * hi) >> 20;
		x[i] = ((lo >> (4 * i)) & 0x0F) + ((hi - 17 * a) << 4);
		hi = a;
	}

	/*
	 * We use the final value for 'hi' as-is; if it is off-range, this
	 * will be detected in the final test.
	 */
	x[7] = (lo >> 28) + (hi << 4);

	/*
	 * Decoding is correct if all decoded values are less than 257.
	 */
	a = (uint32_t)-1;
	for (i = 0; i < 8; i ++) {
		a &= (uint32_t)(x[i] - 257);
	}
	return a >> 31;
}

#if BAT_AVX2

/*
 * Eight decode8() in parallel. 65 bytes are read from buf[], and decoded,
 * to yield 64 values in x[].
 * Return value contains the eight success bits (each as a 32-bit word
 * with value -1 = success, 0 = failure).
 */
TARGET_AVX2
static inline __m256i
decode8x8(uint16_t *x, const uint8_t *buf)
{
	/*
	 * Useful facts:
	 *    256 = 1 mod 17
	 *    65536 = 1 mod 17
	 *    2^32 = 1 mod 17
	 *    for any x in 0..69631, (61681 * x) >> 20 = floor(x / 17)
	 *    17 * 4042322161 = 1 mod 2^32
	 */
	__m256i xa, xb, xm4, xm8, xm16, x257, x61681, x4042322161;
	__m256i xo[8];
	__m256i xi0, xi1, xlo, xhi, xtt;
	__m256i xr0, xr1, xr2, xr3, xs0, xs1, xs2, xs3;
	int i;

	xm4 = _mm256_set1_epi32(0x0F);
	xm8 = _mm256_set1_epi32(0xFF);
	xm16 = _mm256_set1_epi32(0xFFFF);
	x257 = _mm256_set1_epi32(257);
	x61681 = _mm256_set1_epi32(61681);
	x4042322161 = _mm256_set1_epi32(-252645135);

	/* Load low and high parts of each 64-bit word. */
	xi0 = _mm256_loadu_si256((void *)(buf +  0));
	xi1 = _mm256_loadu_si256((void *)(buf + 32));
	xi0 = _mm256_shuffle_epi32(xi0, 0xD8);
	xi1 = _mm256_shuffle_epi32(xi1, 0xD8);
	xi0 = _mm256_permute4x64_epi64(xi0, 0xD8);
	xi1 = _mm256_permute4x64_epi64(xi1, 0xD8);
	xlo = _mm256_permute2x128_si256(xi0, xi1, 0x20);
	xhi = _mm256_permute2x128_si256(xi0, xi1, 0x31);

	/* Load the eight high bits. */
	xtt = _mm256_and_si256(
		_mm256_set1_epi32(1),
		_mm256_srlv_epi32(
			_mm256_set1_epi32(buf[64]),
			_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7)));

	/*
	 * First iteration: we need to compute the quotient and
	 * remainder of (hi + (tt << 32)) divided by 17.
	 *
	 * We first compute the remainder.
	 */
	/* a <- (hi & 0xFFFF) + (hi >> 16) + tt */
	xa = _mm256_add_epi32(
		_mm256_and_si256(xhi, xm16),
		_mm256_add_epi32(_mm256_srli_epi32(xhi, 16), xtt));
	/* a <- (a & 0xFF) + (a >> 8) */
	xa = _mm256_add_epi32(
		_mm256_and_si256(xa, xm8),
		_mm256_srli_epi32(xa, 8));
	/* a <- a - 17*((61681 * a) >> 20) */
	xb = _mm256_srli_epi32(_mm256_mulhi_epu16(xa, x61681), 4);
	xa = _mm256_sub_epi32(xa,
		_mm256_add_epi32(xb, _mm256_slli_epi32(xb, 4)));
	/* x[0] <- (lo & 0x0F) + (a << 4) */
	xo[0] = _mm256_add_epi32(
		_mm256_and_si256(xlo, xm4),
		_mm256_slli_epi32(xa, 4));

	/*
	 * Subtract the remainder, then do the exact division by multiplying
	 * by the inverse of 17 modulo 2^32. We do all computations modulo
	 * 2^32, since we know that the result will be lower than 2^32 (we
	 * don't have to care about the value of tt here, it is implictly
	 * taken into account in the remainder).
	 */
	/* hi <- (hi - a) * 4042322161 */
	xhi = _mm256_mullo_epi32(_mm256_sub_epi32(xhi, xa), x4042322161);

	/*
	 * Next four iterations are similar, except that we work on a
	 * 32-bit value now (no tt).
	 */
	for (i = 1; i <= 4; i ++) {
		/* a <- (hi & 0xFFFF) + (hi >> 16) */
		xa = _mm256_add_epi32(
			_mm256_and_si256(xhi, xm16),
			_mm256_srli_epi32(xhi, 16));
		/* a <- (a & 0xFF) + (a >> 8) */
		xa = _mm256_add_epi32(
			_mm256_and_si256(xa, xm8),
			_mm256_srli_epi32(xa, 8));
		/* a <- a - 17*((61681 * a) >> 20) */
		xb = _mm256_srli_epi32(_mm256_mulhi_epu16(xa, x61681), 4);
		xa = _mm256_sub_epi32(xa,
			_mm256_add_epi32(xb, _mm256_slli_epi32(xb, 4)));
		/* lo <- lo >> 4 */
		xlo = _mm256_srli_epi32(xlo, 4);
		/* x[i] <- (lo & 0x0F) + (a << 4) */
		xo[i] = _mm256_add_epi32(
			_mm256_and_si256(xlo, xm4),
			_mm256_slli_epi32(xa, 4));
		/* hi <- (hi - a) * 4042322161 */
		xhi = _mm256_mullo_epi32(
			_mm256_sub_epi32(xhi, xa), x4042322161);
	}

	/*
	 * At that point, max value for hi is 6049. The next two iterations
	 * can be done with simpler computations.
	 */
	for (i = 5; i <= 6; i ++) {
		/* a <- (61681 * hi) >> 20 */
		xa = _mm256_srli_epi32(_mm256_mulhi_epu16(xhi, x61681), 4);
		/* lo <- lo >> 4 */
		xlo = _mm256_srli_epi32(xlo, 4);
		/* x[i] <- (lo & 0x0F) + ((hi - 17 * a) << 4) */
		xo[i] = _mm256_add_epi32(
			_mm256_and_si256(xlo, xm4),
			_mm256_slli_epi32(
				_mm256_sub_epi32(xhi,
					_mm256_add_epi32(
						xa,
						_mm256_slli_epi32(xa, 4))),
				4));
		/* hi <- a */
		xhi = xa;
	}

	/*
	 * We use the final value for 'hi' as-is; if it is off-range, this
	 * will be detected in the final test.
	 */
	/* lo <- lo >> 4 */
	xlo = _mm256_srli_epi32(xlo, 4);
	/* x[7] <- lo + (hi << 4) */
	xo[7] = _mm256_add_epi32(xlo, _mm256_slli_epi32(xhi, 4));

	/*
	 * Decoding is correct if all decoded values are less than 257.
	 */
	xa = _mm256_and_si256(
		_mm256_and_si256(
			_mm256_and_si256(
				_mm256_cmpgt_epi32(x257, xo[0]),
				_mm256_cmpgt_epi32(x257, xo[1])),
			_mm256_and_si256(
				_mm256_cmpgt_epi32(x257, xo[2]),
				_mm256_cmpgt_epi32(x257, xo[3]))),
		_mm256_and_si256(
			_mm256_and_si256(
				_mm256_cmpgt_epi32(x257, xo[4]),
				_mm256_cmpgt_epi32(x257, xo[5])),
			_mm256_and_si256(
				_mm256_cmpgt_epi32(x257, xo[6]),
				_mm256_cmpgt_epi32(x257, xo[7]))));

	/*
	 * Encode values. Registers currently contain (order low -> high
	 * in each register, 16-bit elements):
	 *   xo[0]   x00 0 x08 0 x16 0 x24 0 x32 0 x40 0 x48 0 x56 0
	 *   xo[1]   x01 0 x09 0 x17 0 x25 0 x33 0 x41 0 x49 0 x57 0
	 *   xo[2]   x02 0 x10 0 x18 0 x26 0 x34 0 x42 0 x50 0 x58 0
	 *   xo[3]   x03 0 x11 0 x19 0 x27 0 x35 0 x43 0 x51 0 x59 0
	 *   xo[4]   x04 0 x12 0 x20 0 x28 0 x36 0 x44 0 x52 0 x60 0
	 *   xo[5]   x05 0 x13 0 x21 0 x29 0 x37 0 x45 0 x53 0 x61 0
	 *   xo[6]   x06 0 x14 0 x22 0 x30 0 x38 0 x46 0 x54 0 x62 0
	 *   xo[7]   x07 0 x15 0 x23 0 x31 0 x39 0 x47 0 x55 0 x63 0
	 * All the return values fit on 16 bits; we ignore the high half
	 * of each 32-bit element.
	 *
	 * We first pack values into four registers.
	 */
	xr0 = _mm256_or_si256(xo[0], _mm256_slli_epi32(xo[1], 16));
	xr1 = _mm256_or_si256(xo[2], _mm256_slli_epi32(xo[3], 16));
	xr2 = _mm256_or_si256(xo[4], _mm256_slli_epi32(xo[5], 16));
	xr3 = _mm256_or_si256(xo[6], _mm256_slli_epi32(xo[7], 16));

	/*
	 * xr0: x00 x01 x08 x09 x16 x17 x24 x25 x32 x33 x40 x41 x48 x49 x56 x57
	 * xr1: x02 x03 x10 x11 x18 x19 x26 x27 x34 x35 x42 x43 x50 x51 x58 x59
	 * xr2: x04 x05 x12 x13 x20 x21 x28 x29 x36 x37 x44 x45 x52 x53 x60 x61
	 * xr3: x06 x07 x14 x15 x22 x23 x30 x31 x38 x39 x46 x47 x54 x55 x62 x63
	 */

	xs0 = _mm256_unpacklo_epi32(xr0, xr1);
	xs1 = _mm256_unpackhi_epi32(xr0, xr1);
	xs2 = _mm256_unpacklo_epi32(xr2, xr3);
	xs3 = _mm256_unpackhi_epi32(xr2, xr3);

	/*
	 * xs0: x00 x01 x02 x03 x08 x09 x10 x11 x32 x33 x34 x35 x40 x41 x42 x43
	 * xs1: x16 x17 x18 x19 x24 x25 x26 x27 x48 x19 x50 x51 x56 x57 x58 x59
	 * xs2: x04 x05 x06 x07 x12 x13 x14 x15 x36 x37 x38 x39 x44 x45 x46 x47
	 * xs3: x20 x21 x22 x23 x28 x29 x30 x31 x52 x53 x54 x55 x60 x61 x62 x63
	 */

	xr0 = _mm256_unpacklo_epi64(xs0, xs2);
	xr1 = _mm256_unpackhi_epi64(xs0, xs2);
	xr2 = _mm256_unpacklo_epi64(xs1, xs3);
	xr3 = _mm256_unpackhi_epi64(xs1, xs3);

	/*
	 * xr0: x00 x01 x02 x03 x04 x05 x06 x07 x32 x33 x34 x35 x36 x37 x38 x39
	 * xr1: x08 x09 x10 x11 x12 x13 x14 x15 x40 x41 x42 x43 x44 x45 x46 x47
	 * xr2: x16 x17 x18 x19 x20 x21 x22 x23 x48 x49 x50 x51 x52 x53 x54 x55
	 * xr3: x24 x25 x26 x27 x28 x29 x30 x31 x56 x57 x58 x59 x60 x61 x62 x63
	 */

	xs0 = _mm256_permute2x128_si256(xr0, xr1, 0x20);
	xs1 = _mm256_permute2x128_si256(xr2, xr3, 0x20);
	xs2 = _mm256_permute2x128_si256(xr0, xr1, 0x31);
	xs3 = _mm256_permute2x128_si256(xr2, xr3, 0x31);

	_mm256_storeu_si256((void *)(x +  0), xs0);
	_mm256_storeu_si256((void *)(x + 16), xs1);
	_mm256_storeu_si256((void *)(x + 32), xs2);
	_mm256_storeu_si256((void *)(x + 48), xs3);

	return xa;
}

#endif

/* see inner.h */
size_t
bat_encode_257(void *out, size_t max_out_len,
	const uint16_t *x, unsigned logn)
{
	size_t out_len, u, v, n;
	uint8_t *buf;
	uint32_t xb;
	int j;

	switch (logn) {
	case 1:
		out_len = 3;
		break;
	case 2:
		out_len = 5;
		break;
	default:
		out_len = (((1u << (logn - 3)) + 7) >> 3) + (1u << logn);
		break;
	}
	if (out == NULL) {
		return out_len;
	}
	if (max_out_len < out_len) {
		return 0;
	}

	buf = out;
	switch (logn) {
	case 1:
		buf[0] = (uint8_t)x[0];
		buf[1] = (uint8_t)((x[0] >> 8) | ((unsigned)x[1] << 2));
		buf[2] = (uint8_t)(x[1] >> 6);
		return 3;
	case 2:
		buf[0] = (uint8_t)x[0];
		buf[1] = (uint8_t)((x[0] >> 8) | ((unsigned)x[1] << 2));
		buf[2] = (uint8_t)((x[1] >> 6) | ((unsigned)x[2] << 4));
		buf[3] = (uint8_t)((x[2] >> 4) | ((unsigned)x[3] << 6));
		buf[4] = (uint8_t)(x[3] >> 2);
		return 5;
	case 3:
	case 4:
	case 5:
		n = (size_t)1 << logn;
		v = 0;
		xb = 0;
		j = 0;
		for (u = 0; u < n; u += 8) {
			uint32_t lo, hi;

			xb |= encode8(&lo, &hi, x + u) << j;
			enc32le(buf + v, lo);
			enc32le(buf + v + 4, hi);
			v += 8;
			if (++ j == 8) {
				buf[v ++] = xb;
				xb = 0;
				j = 0;
			}
		}
		if (j != 0) {
			buf[v ++] = xb;
		}
		return v;
	default:
		n = (size_t)1 << logn;
		v = 0;
#if BAT_AVX2
		for (u = 0; u < n; u += 64) {
			buf[v + 64] = encode8x8(buf + v, x + u);
			v += 65;
		}
#else
		for (u = 0; u < n; u += 64) {
			xb = 0;
			for (j = 0; j < 8; j ++) {
				uint32_t lo, hi;

				xb |= encode8(&lo, &hi, x + u + (j << 3)) << j;
				enc32le(buf + v + (j << 3), lo);
				enc32le(buf + v + (j << 3) + 4, hi);
			}
			buf[v + 64] = xb;
			v += 65;
		}
#endif
		return v;
	}
}

/* see inner.h */
size_t
bat_decode_257(uint16_t *x, unsigned logn,
	const void *in, size_t max_in_len)
{
	const uint8_t *buf;
	uint32_t r, xb;
	size_t in_len, u, n;
#if BAT_AVX2
	__m256i xr;
#endif

	buf = in;
	switch (logn) {
	case 1:
		if (max_in_len < 3) {
			return 0;
		}
		x[0] = buf[0] | ((buf[1] & 0x03) << 8);
		x[1] = (buf[1] >> 2) | (buf[2] << 6);
		in_len = 3;
		r = (((uint32_t)x[0] - 257) & ((uint32_t)x[1] - 257)) >> 31;
		break;
	case 2:
		if (max_in_len < 5) {
			return 0;
		}
		x[0] = buf[0] | ((buf[1] & 0x03) << 8);
		x[1] = (buf[1] >> 2) | ((buf[2] & 0x0F) << 6);
		x[2] = (buf[2] >> 4) | ((buf[3] & 0x3F) << 4);
		x[3] = (buf[3] >> 6) | (buf[4] << 2);
		in_len = 5;
		r = (((uint32_t)x[0] - 257)
			& ((uint32_t)x[1] - 257)
			& ((uint32_t)x[2] - 257)
			& ((uint32_t)x[3] - 257)) >> 31;
		break;
	case 3:
	case 4:
	case 5:
		n = (size_t)1 << logn;
		in_len = n + 1;
		if (max_in_len < in_len) {
			return 0;
		}
		xb = buf[n];

		/* Ignored bits in xb must be zero. */
		r = (uint32_t)((xb >> (1u << (logn - 3))) - 1) >> 31;

		/* Decode 65-bit blocks. */
		for (u = 0; u < n; u += 8) {
			r &= decode8(x + u,
				dec32le(buf + u),
				dec32le(buf + u + 4),
				(xb >> ((unsigned)u >> 3)) & 1);
		}
		break;

	default:
		n = (size_t)1 << logn;
		in_len = n + ((size_t)1 << (logn - 6));
		if (max_in_len < in_len) {
			return 0;
		}
#if BAT_AVX2
		xr = _mm256_set1_epi32(-1);
		for (u = 0; u < n; u += 64) {
			xr = _mm256_and_si256(xr, decode8x8(x + u, buf));
			buf += 65;
		}

		xr = _mm256_and_si256(xr, _mm256_bsrli_epi128(xr, 4));
		xr = _mm256_and_si256(xr, _mm256_bsrli_epi128(xr, 8));
		r = _mm_cvtsi128_si32(
			_mm_and_si128(
				_mm256_castsi256_si128(xr),
				_mm256_extracti128_si256(xr, 1))) & 1;
#else
		r = 1;
		for (u = 0; u < n; u += 64) {
			size_t v;

			xb = buf[64];
			for (v = 0; v < 64; v += 8) {
				r &= decode8(x + u + v,
					dec32le(buf + v),
					dec32le(buf + v + 4),
					(xb >> ((unsigned)v >> 3)) & 1);
			}
			buf += 65;
		}
#endif
		break;
	}

	return in_len & -(size_t)r;
}

/*
 * Ciphertext encoding: this is similar to public keys, except that
 * values are in the -64..+64 range. We first add 64 to get a positive
 * value (0..128). Each value will be split into its low 3 bits and high
 * 5 bits; the high bits are in 0..16 only. Values are grouped into
 * chunks of eight, that encode into 7 bytes + 1 bit.
 *
 * For small degrees:
 *
 *  - If degree is 2, then encoding is over two bytes, 8 bits per value.
 *  - If degree is 4, then encoding is over four bytes, 8 bits per value.
 *  - If degree is 8, 16 or 32, then one, two or four 57-bit blocks are
 *    used; the ignored bits in the last byte are set to zero.
 */

/*
 * Encode eight ciphertext values x[0]..x[7] into 57 bits (bits 0..23 are
 * returned in *lo, bits 24..55 in *hi, bit 56 is the function return value).
 */
static inline uint32_t
encode_ct_8(uint32_t *lo, uint32_t *hi, const int8_t *x)
{
	int i;
	uint32_t wl, wh, r, t;

	wl = 0;
	for (i = 0; i < 8; i ++) {
		wl |= (uint32_t)((uint8_t)x[i] & 0x07) << (3 * i);
	}
	*lo = wl;
	wh = 0;
	for (i = 7; i > 0; i --) {
		wh = (wh * 17) + (uint32_t)((x[i] + 64) >> 3);
	}

	/*
	 * Final iteration is:
	 *   *hi = wh*17 + (x[0] >> 3)
	 * which can be written as:
	 *   *hi = wh*16 + wh + (x[0] >> 3)
	 * There will be a carry into the top bit (bit 64 of the final
	 * output) if and only if one of the following holds:
	 *  - wh >= 2^28
	 *  - adding wh + (x[0] >> 3) triggers a carry
	 * Note that wh + (x[0] >> 3) <= 410338688 < 2^29.
	 */

	r = wh >> 28;
	t = wh + ((x[0] + 64) >> 3);
	wh = (wh << 4) + t;
	r |= ((uint32_t)(wh - t) & ~wh) >> 31;
	*hi = wh;
	return r;
}

#if BAT_AVX2

/*
 * Eight encode_ct_8() in parallel. 64 input values (x[0..63]) are encoded
 * into 56 consecutive bytes (buf[0..55]) and the eight extra bits are
 * returned as a packed byte value.
 */
TARGET_AVX2
static inline uint32_t
encode_ct_8x8(uint8_t *buf, const int8_t *x)
{
	__m256i xm3, xm8, xm16, xm32, x64, x83521;
	__m256i xi0, xi1, xa0, xa1, xt0, xt1;

	xm3 = _mm256_set1_epi8(0x07);
	xm8 = _mm256_set1_epi16(0x00FF);
	xm16 = _mm256_set1_epi32(0xFFFF);
	xm32 = _mm256_set1_epi64x(0xFFFFFFFF);
	x64 = _mm256_set1_epi8(64);
	x83521 = _mm256_set1_epi64x(83521);

	/*
	 * Load the 64 values to encode and normalize into 0..128.
	 */
	xi0 = _mm256_loadu_si256((const void *)(x +  0));
	xi1 = _mm256_loadu_si256((const void *)(x + 32));
	xi0 = _mm256_add_epi8(xi0, x64);
	xi1 = _mm256_add_epi8(xi1, x64);

	/*
	 * Low parts: we extract the low 3 bits of each value, and
	 * repack them.
	 */
	xa0 = _mm256_and_si256(xi0, xm3);
	xa1 = _mm256_and_si256(xi1, xm3);

	xa0 = _mm256_or_si256(
		_mm256_and_si256(xa0, xm8),
		_mm256_srli_epi16(xa0, 5));
	xa1 = _mm256_or_si256(
		_mm256_and_si256(xa1, xm8),
		_mm256_srli_epi16(xa1, 5));

	xa0 = _mm256_or_si256(
		_mm256_and_si256(xa0, xm16),
		_mm256_srli_epi32(xa0, 10));
	xa1 = _mm256_or_si256(
		_mm256_and_si256(xa1, xm16),
		_mm256_srli_epi32(xa1, 10));

	xa0 = _mm256_or_si256(
		_mm256_and_si256(xa0, xm32),
		_mm256_srli_epi64(xa0, 20));
	xa1 = _mm256_or_si256(
		_mm256_and_si256(xa1, xm32),
		_mm256_srli_epi64(xa1, 20));

	/*
	 * Low parts are now in xa0 and xa1 (24-bit values, each in a
	 * 32-bit element):
	 *
	 *   xa0:  lo0 0 lo1 0 lo2 0 lo3 0
	 *   xa1:  lo4 0 lo5 0 lo6 0 lo7 0
	 */

	/*
	 * Extract high parts from each input value.
	 */
	xi0 = _mm256_srli_epi16(_mm256_andnot_si256(xm3, xi0), 3);
	xi1 = _mm256_srli_epi16(_mm256_andnot_si256(xm3, xi1), 3);

	/*
	 * Combine high parts by pairs (one 16-bit element per pair).
	 */
	xi0 = _mm256_add_epi16(
		_mm256_srli_epi16(_mm256_andnot_si256(xm8, xi0), 4),
		_mm256_add_epi16(
			_mm256_and_si256(xi0, xm8),
			_mm256_srli_epi16(xi0, 8)));
	xi1 = _mm256_add_epi16(
		_mm256_srli_epi16(_mm256_andnot_si256(xm8, xi1), 4),
		_mm256_add_epi16(
			_mm256_and_si256(xi1, xm8),
			_mm256_srli_epi16(xi1, 8)));

	/*
	 * Combine pairs into quartets (one 32-bit element per quartet).
	 * The high pair is of each quartet is multiplied by 289, then
	 * added to the low pair. Since multiplications have relatively
	 * high latency, we use two "manual" multiplications by 17 instead.
	 */
	xt0 = _mm256_srli_epi32(xi0, 16);
	xt1 = _mm256_srli_epi32(xi1, 16);
	xt0 = _mm256_add_epi32(xt0, _mm256_slli_epi32(xt0, 4));
	xt1 = _mm256_add_epi32(xt1, _mm256_slli_epi32(xt1, 4));
	xt0 = _mm256_add_epi32(xt0, _mm256_slli_epi32(xt0, 4));
	xt1 = _mm256_add_epi32(xt1, _mm256_slli_epi32(xt1, 4));
	xi0 = _mm256_add_epi32(_mm256_and_si256(xi0, xm16), xt0);
	xi1 = _mm256_add_epi32(_mm256_and_si256(xi1, xm16), xt1);

	/*
	 * Combine quartets together. This time, mul_epu32() is more
	 * efficient than four manual multiplications by 17. We get
	 * each result (high parts, 33 bits) into a 64-bit element.
	 */
	xi0 = _mm256_add_epi64(
		_mm256_and_si256(xi0, xm32),
		_mm256_mul_epu32(_mm256_srli_epi64(xi0, 32), x83521));
	xi1 = _mm256_add_epi64(
		_mm256_and_si256(xi1, xm32),
		_mm256_mul_epu32(_mm256_srli_epi64(xi1, 32), x83521));

	/*
	 * High parts (hi*, 32-bit) and extra bits (t*, 1-bit) are now in
	 * xi0 and xi1:
	 *
	 *   xi0: hi0 t0 hi1 t1 hi2 t2 hi3 t3
	 *   xi1: hi4 t4 hi5 t5 hi6 t6 hi7 t7
	 */

	/*
	 * Combine the low (24-bit) and high (32-bit) parts (without the
	 * t* bits) into 56 bytes. We first make the eight 7-byte chunks
	 * in 64-bit elements.
	 */
	xa0 = _mm256_or_si256(xa0, _mm256_slli_epi64(xi0, 24));
	xa1 = _mm256_or_si256(xa1, _mm256_slli_epi64(xi1, 24));

	/*
	 * We now need to make four output values, of 4, 4, 16 and 32
	 * bytes, respectively, so that we may do non-overlapping writes
	 * (to avoid aliasing issues) and avoid writing beyond the
	 * output buffer. The first eight bytes need to be written as
	 * two chunks of 4 bytes because we might be compiling in 32-bit
	 * mode, in which case intrinsics to extract the low 64 bits as
	 * an integer may not be available.
	 *
	 */

	/* Output 0..3 */
	enc32le(buf, (uint32_t)_mm256_cvtsi256_si32(xa0));

	/* Output 4..7 */
	enc32le(buf + 4, (uint32_t)_mm_cvtsi128_si32(
		_mm_shuffle_epi8(_mm256_castsi256_si128(xa0),
			_mm_setr_epi8(4, 5, 6, 8, -1, -1, -1, -1,
				-1, -1, -1, -1, -1, -1, -1, -1))));

	/* Output 8..23 */
	xt0 = _mm256_bsrli_epi128(xa0, 9);
	xt1 = _mm256_shuffle_epi8(xa0,
		_mm256_setr_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1,
			0, 1, 2, 3, 4, 5, 6, 8, 9, 10));
	_mm_storeu_si128((void *)(buf + 8),
		_mm_blend_epi16(
			_mm256_castsi256_si128(xt0),
			_mm256_extracti128_si256(xt1, 1),
			0xF8));

	/* Output 24..55 */
	/* We want 4 bytes from the high lane of xa0 (into the low
	   lane of the result), 14 bytes from the low lane of xa1
	   (12 into the low lane of the result, 2 into the high lane)
	   and 14 bytes from the high lane of xa1 (all into the high
	   lane of the result). We only need one cross-lane move, which
	   we do first (into xt0). */
	xt0 = _mm256_permute2x128_si256(xa0, xa1, 0x21);
	xt1 = _mm256_or_si256(
		_mm256_or_si256(
			_mm256_shuffle_epi8(xt0, _mm256_setr_epi8(
				11, 12, 13, 14, -1, -1, -1, -1,
				-1, -1, -1, -1, -1, -1, -1, -1,
				-1, -1, -1, -1, -1, -1, -1, -1,
				-1, -1, -1, -1, -1, -1, -1, -1)),
			_mm256_shuffle_epi8(xa1, _mm256_setr_epi8(
				-1, -1, -1, -1,
				0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12,
				-1, -1, -1, -1, -1, -1, -1, -1,
				-1, -1, -1, -1, -1, -1, -1, -1))),
		_mm256_or_si256(
			_mm256_shuffle_epi8(xt0, _mm256_setr_epi8(
				-1, -1, -1, -1, -1, -1, -1, -1,
				-1, -1, -1, -1, -1, -1, -1, -1,
				13, 14, -1, -1, -1, -1, -1, -1,
				-1, -1, -1, -1, -1, -1, -1, -1)),
			_mm256_shuffle_epi8(xa1, _mm256_setr_epi8(
				-1, -1, -1, -1, -1, -1, -1, -1,
				-1, -1, -1, -1, -1, -1, -1, -1,
				-1, -1, 0, 1, 2, 3, 4, 5,
				6, 8, 9, 10, 11, 12, 13, 14))));
	_mm256_storeu_si256((void *)(buf + 24), xt1);

	/*
	 * Extract the extra bits to make the extra byte.
	 */
	xi0 = _mm256_or_si256(
		_mm256_srli_epi64(xi0, 32),
		_mm256_slli_epi64(_mm256_bsrli_epi128(xi0, 12), 1));
	xi1 = _mm256_or_si256(
		_mm256_srli_epi64(xi1, 32),
		_mm256_slli_epi64(_mm256_bsrli_epi128(xi1, 12), 1));
	xi0 = _mm256_or_si256(xi0, _mm256_slli_epi64(xi1, 4));
	xi0 = _mm256_permute4x64_epi64(xi0, 0xD8);
	xi0 = _mm256_or_si256(xi0, _mm256_bsrli_epi128(xi0, 1));
	xi0 = _mm256_or_si256(xi0, _mm256_srli_epi64(xi0, 54));
	return _mm256_cvtsi256_si32(xi0);
}

#endif

/*
 * Decode a 57-bit word (lo is bits 0..23, hi is bits 24..55, tt is bit 56)
 * into eight values. Returned value is 1 on success, 0 on error.
 */
static uint32_t
decode_ct_8(int8_t *x, uint32_t lo, uint32_t hi, uint32_t tt)
{
	/*
	 * Useful facts:
	 *    256 = 1 mod 17
	 *    65536 = 1 mod 17
	 *    2^32 = 1 mod 17
	 *    for any x in 0..69631, (61681 * x) >> 20 = floor(x / 17)
	 *    17 * 4042322161 = 1 mod 2^32
	 */
	uint32_t a;
	int i;

	/*
	 * First iteration: we need to compute the quotient and
	 * remainder of (hi + (tt << 32)) divided by 17.
	 *
	 * We first compute the remainder.
	 */
	a = (hi & 0xFFFF) + (hi >> 16) + tt;   /* 0 <= a <= 131071 */
	a = (a & 0xFF) + (a >> 8);             /* 0 <= a <= 766 */
	a -= 17 * ((61681 * a) >> 20);         /* 0 <= a <= 16 */
	x[0] = (lo & 0x07) + (a << 3) - 64;

	/*
	 * Subtract the remainder, then do the exact division by multiplying
	 * by the inverse of 17 modulo 2^32. We do all computations modulo
	 * 2^32, since we know that the result will be lower than 2^32 (we
	 * don't have to care about the value of tt here, it is implictly
	 * taken into account in the remainder).
	 */
	hi = (hi - a) * 4042322161u;

	/*
	 * Next four iterations are similar, except that we work on a
	 * 32-bit value now (no tt).
	 */
	for (i = 1; i <= 4; i ++) {
		a = (hi & 0xFFFF) + (hi >> 16);   /* 0 <= a <= 131070 */
		a = (a & 0xFF) + (a >> 8);        /* 0 <= a <= 766 */
		a -= 17 * ((61681 * a) >> 20);    /* 0 <= a <= 16 */
		x[i] = ((lo >> (3 * i)) & 0x07) + (a << 3) - 64;
		hi = (hi - a) * 4042322161u;
	}

	/*
	 * At that point, max value for hi is 6049. The next two iterations
	 * can be done with simpler computations.
	 */
	for (i = 5; i <= 6; i ++) {
		a = (61681 * hi) >> 20;
		x[i] = ((lo >> (3 * i)) & 0x07) + ((hi - 17 * a) << 3) - 64;
		hi = a;
	}

	/*
	 * We use the final value for 'hi' as-is; if it is off-range, this
	 * will be detected in the final test.
	 */
	x[7] = (lo >> 21) + (hi << 3) - 64;

	/*
	 * Decoding is correct if all decoded values are in -64..+64.
	 * Since the values were computed by subtracting 64 from a nonnegative
	 * integer, they cannot be lower than -64, but they could be higher
	 * than +64.
	 */
	a = (uint32_t)-1;
	for (i = 0; i < 8; i ++) {
		a &= (uint32_t)(x[i] - 65);
	}
	return a >> 31;
}

#if BAT_AVX2

/*
 * Eight decode_ct_8() in parallel. 57 bytes are read from buf[] and
 * decoded into 64 values into x[].
 * Return value contains the eight success bits as 32-bit elements
 * (with value -1 = success, 0 = failure).
 */
TARGET_AVX2
static inline __m256i
decode_ct_8x8(int8_t *x, const uint8_t *buf)
{
	/*
	 * Useful facts:
	 *    256 = 1 mod 17
	 *    65536 = 1 mod 17
	 *    2^32 = 1 mod 17
	 *    for any x in 0..69631, (61681 * x) >> 20 = floor(x / 17)
	 *    17 * 4042322161 = 1 mod 2^32
	 */
	__m256i xa, xb, xm3, xm8, xm16, x64, x129, x61681, x4042322161;
	__m256i xo[8];
	__m256i xi0, xi1, xlo, xhi, xtt;
	__m256i xr0, xr1, xs0, xs1;
	int i;

	xm3 = _mm256_set1_epi32(0x07);
	xm8 = _mm256_set1_epi32(0xFF);
	xm16 = _mm256_set1_epi32(0xFFFF);
	x64 = _mm256_set1_epi8(64);
	x129 = _mm256_set1_epi32(129);
	x61681 = _mm256_set1_epi32(61681);
	x4042322161 = _mm256_set1_epi32(-252645135);

	/* Load the 56 bytes. */
	xi0 = _mm256_loadu_si256((void *)(buf + 0));
	xi1 = _mm256_loadu_si256((void *)(buf + 24));

	/* Move bytes around so that we get the low parts (24 bits each)
	   into xlo (in 32-bit elements) and the high parts (32 bits each)
	   into xhi. */

	/* xi0 <- bytes 0..15 | 16..31
	   xi1 <- bytes 24..39 | 40..55
	   xa  <- bytes 16..31 | 24..39 */
	xa = _mm256_permute2x128_si256(xi0, xi1, 0x21);

	/* xlo <- bytes:
	    0 1 2 . 7 8 9 . 14 15 16 . 21 22 23 .
	    28 29 30 . 35 36 37 . 42 43 44 . 49 50 51 . */
	xlo = _mm256_or_si256(
		_mm256_or_si256(
			_mm256_shuffle_epi8(xi0, _mm256_setr_epi8(
				0, 1, 2, -1, 7, 8, 9, -1,
				14, 15, -1, -1, -1, -1, -1, -1,
				12, 13, 14, -1, -1, -1, -1, -1,
				-1, -1, -1, -1, -1, -1, -1, -1)),
			_mm256_shuffle_epi8(xi1, _mm256_setr_epi8(
				-1, -1, -1, -1, -1, -1, -1, -1,
				-1, -1, -1, -1, -1, -1, -1, -1,
				-1, -1, -1, -1, -1, -1, -1, -1,
				2, 3, 4, -1, 9, 10, 11, -1))),
		_mm256_shuffle_epi8(xa, _mm256_setr_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, 0, -1, 5, 6, 7, -1,
			-1, -1, -1, -1, 11, 12, 13, -1,
			-1, -1, -1, -1, -1, -1, -1, -1)));

	/* xhi <- bytes:
	    3 4 5 6 10 11 12 13 17 18 19 20 24 25 26 27
	    31 32 33 34 38 39 40 41 45 46 47 48 52 53 54 55 */
	xhi = _mm256_or_si256(
		_mm256_or_si256(
			_mm256_shuffle_epi8(xi0, _mm256_setr_epi8(
				3, 4, 5, 6, 10, 11, 12, 13,
				-1, -1, -1, -1, -1, -1, -1, -1,
				15, -1, -1, -1, -1, -1, -1, -1,
				-1, -1, -1, -1, -1, -1, -1, -1)),
			_mm256_shuffle_epi8(xi1, _mm256_setr_epi8(
				-1, -1, -1, -1, -1, -1, -1, -1,
				-1, -1, -1, -1, 0, 1, 2, 3,
				-1, -1, -1, -1, -1, -1, 0, 1,
				5, 6, 7, 8, 12, 13, 14, 15))),
		_mm256_shuffle_epi8(xa, _mm256_setr_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1,
			1, 2, 3, 4, -1, -1, -1, -1,
			-1, 8, 9, 10, 14, 15, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1)));

	/* Extract the high bits into xtt. */
	xtt = _mm256_and_si256(
		_mm256_set1_epi32(1),
		_mm256_srlv_epi32(
			_mm256_set1_epi32(buf[56]),
			_mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7)));

	/*
	 * First iteration: we need to compute the quotient and
	 * remainder of (hi + (tt << 32)) divided by 17.
	 *
	 * We first compute the remainder.
	 */
	/* a <- (hi & 0xFFFF) + (hi >> 16) + tt */
	xa = _mm256_add_epi32(
		_mm256_and_si256(xhi, xm16),
		_mm256_add_epi32(_mm256_srli_epi32(xhi, 16), xtt));
	/* a <- (a & 0xFF) + (a >> 8) */
	xa = _mm256_add_epi32(
		_mm256_and_si256(xa, xm8),
		_mm256_srli_epi32(xa, 8));
	/* a <- a - 17*((61681 * a) >> 20) */
	xb = _mm256_srli_epi32(_mm256_mulhi_epu16(xa, x61681), 4);
	xa = _mm256_sub_epi32(xa,
		_mm256_add_epi32(xb, _mm256_slli_epi32(xb, 4)));
	/* x[0] <- (lo & 0x07) + (a << 3) */
	xo[0] = _mm256_add_epi32(
		_mm256_and_si256(xlo, xm3),
		_mm256_slli_epi32(xa, 3));

	/*
	 * Subtract the remainder, then do the exact division by multiplying
	 * by the inverse of 17 modulo 2^32. We do all computations modulo
	 * 2^32, since we know that the result will be lower than 2^32 (we
	 * don't have to care about the value of tt here, it is implictly
	 * taken into account in the remainder).
	 */
	/* hi <- (hi - a) * 4042322161 */
	xhi = _mm256_mullo_epi32(_mm256_sub_epi32(xhi, xa), x4042322161);

	/*
	 * Next four iterations are similar, except that we work on a
	 * 32-bit value now (no tt).
	 */
	for (i = 1; i <= 4; i ++) {
		/* a <- (hi & 0xFFFF) + (hi >> 16) */
		xa = _mm256_add_epi32(
			_mm256_and_si256(xhi, xm16),
			_mm256_srli_epi32(xhi, 16));
		/* a <- (a & 0xFF) + (a >> 8) */
		xa = _mm256_add_epi32(
			_mm256_and_si256(xa, xm8),
			_mm256_srli_epi32(xa, 8));
		/* a <- a - 17*((61681 * a) >> 20) */
		xb = _mm256_srli_epi32(_mm256_mulhi_epu16(xa, x61681), 4);
		xa = _mm256_sub_epi32(xa,
			_mm256_add_epi32(xb, _mm256_slli_epi32(xb, 4)));
		/* lo <- lo >> 3 */
		xlo = _mm256_srli_epi32(xlo, 3);
		/* x[i] <- (lo & 0x07) + (a << 3) */
		xo[i] = _mm256_add_epi32(
			_mm256_and_si256(xlo, xm3),
			_mm256_slli_epi32(xa, 3));
		/* hi <- (hi - a) * 4042322161 */
		xhi = _mm256_mullo_epi32(
			_mm256_sub_epi32(xhi, xa), x4042322161);
	}

	/*
	 * At that point, max value for hi is 6049. The next two iterations
	 * can be done with simpler computations.
	 */
	for (i = 5; i <= 6; i ++) {
		/* a <- (61681 * hi) >> 20 */
		xa = _mm256_srli_epi32(_mm256_mulhi_epu16(xhi, x61681), 4);
		/* lo <- lo >> 3 */
		xlo = _mm256_srli_epi32(xlo, 3);
		/* x[i] <- (lo & 0x07) + ((hi - 17 * a) << 3) */
		xo[i] = _mm256_add_epi32(
			_mm256_and_si256(xlo, xm3),
			_mm256_slli_epi32(
				_mm256_sub_epi32(xhi,
					_mm256_add_epi32(
						xa,
						_mm256_slli_epi32(xa, 4))),
				3));
		/* hi <- a */
		xhi = xa;
	}

	/*
	 * We use the final value for 'hi' as-is; if it is off-range, this
	 * will be detected in the final test.
	 */
	/* lo <- lo >> 3 */
	xlo = _mm256_srli_epi32(xlo, 3);
	/* x[7] <- lo + (hi << 3) */
	xo[7] = _mm256_add_epi32(xlo, _mm256_slli_epi32(xhi, 3));

	/*
	 * Decoding is correct if all decoded values are less than 129.
	 */
	xa = _mm256_and_si256(
		_mm256_and_si256(
			_mm256_and_si256(
				_mm256_cmpgt_epi32(x129, xo[0]),
				_mm256_cmpgt_epi32(x129, xo[1])),
			_mm256_and_si256(
				_mm256_cmpgt_epi32(x129, xo[2]),
				_mm256_cmpgt_epi32(x129, xo[3]))),
		_mm256_and_si256(
			_mm256_and_si256(
				_mm256_cmpgt_epi32(x129, xo[4]),
				_mm256_cmpgt_epi32(x129, xo[5])),
			_mm256_and_si256(
				_mm256_cmpgt_epi32(x129, xo[6]),
				_mm256_cmpgt_epi32(x129, xo[7]))));

	/*
	 * Encode values. Registers currently contain (order low -> high
	 * in each register, 16-bit elements):
	 *   xo[0]   x00 0 x08 0 x16 0 x24 0 x32 0 x40 0 x48 0 x56 0
	 *   xo[1]   x01 0 x09 0 x17 0 x25 0 x33 0 x41 0 x49 0 x57 0
	 *   xo[2]   x02 0 x10 0 x18 0 x26 0 x34 0 x42 0 x50 0 x58 0
	 *   xo[3]   x03 0 x11 0 x19 0 x27 0 x35 0 x43 0 x51 0 x59 0
	 *   xo[4]   x04 0 x12 0 x20 0 x28 0 x36 0 x44 0 x52 0 x60 0
	 *   xo[5]   x05 0 x13 0 x21 0 x29 0 x37 0 x45 0 x53 0 x61 0
	 *   xo[6]   x06 0 x14 0 x22 0 x30 0 x38 0 x46 0 x54 0 x62 0
	 *   xo[7]   x07 0 x15 0 x23 0 x31 0 x39 0 x47 0 x55 0 x63 0
	 * All the return values fit on 8 bits; thus, we first pack
	 * values into two registers.
	 */
	xr0 = _mm256_or_si256(
		_mm256_or_si256(
			xo[0],
			_mm256_slli_epi32(xo[1], 8)),
		_mm256_or_si256(
			_mm256_slli_epi32(xo[2], 16),
			_mm256_slli_epi32(xo[3], 24)));
	xr1 = _mm256_or_si256(
		_mm256_or_si256(
			xo[4],
			_mm256_slli_epi32(xo[5], 8)),
		_mm256_or_si256(
			_mm256_slli_epi32(xo[6], 16),
			_mm256_slli_epi32(xo[7], 24)));

	/*
	 * xr0: 0-3  8-11 16-19 24-27 | 32-35 40-43 48-51 56-59
	 * xr1: 4-7 12-15 20-23 38-31 | 36-39 44-47 52-55 60-63
	 */

	xs0 = _mm256_unpacklo_epi32(xr0, xr1);
	xs1 = _mm256_unpackhi_epi32(xr0, xr1);

	/*
	 * xs0:  0-15 | 32-47
	 * xs1: 16-31 | 48-63
	 */

	xr0 = _mm256_permute2x128_si256(xs0, xs1, 0x20);
	xr1 = _mm256_permute2x128_si256(xs0, xs1, 0x31);

	/*
	 * We got values in the 0..128 range, we must normalize to -64..+64.
	 */
	xr0 = _mm256_sub_epi8(xr0, x64);
	xr1 = _mm256_sub_epi8(xr1, x64);

	_mm256_storeu_si256((void *)(x +  0), xr0);
	_mm256_storeu_si256((void *)(x + 32), xr1);

	return xa;
}

#endif

/* see inner.h */
size_t
bat_encode_ciphertext_257(void *out, size_t max_out_len,
	const int8_t *c, unsigned logn)
{
	size_t out_len, u, v, n;
	uint8_t *buf;
	uint32_t xb;
	int j;

	switch (logn) {
	case 1:
		out_len = 2;
		break;
	case 2:
		out_len = 4;
		break;
	default:
		out_len = (((1u << (logn - 3)) + 7) >> 3) + (1u << logn)
			- (1u << (logn - 3));
		break;
	}
	if (out == NULL) {
		return out_len;
	}
	if (max_out_len < out_len) {
		return 0;
	}

	buf = out;
	switch (logn) {
	case 1:
		buf[0] = (uint8_t)(c[0] + 64);
		buf[1] = (uint8_t)(c[1] + 64);
		return 2;
	case 2:
		buf[0] = (uint8_t)(c[0] + 64);
		buf[1] = (uint8_t)(c[1] + 64);
		buf[2] = (uint8_t)(c[2] + 64);
		buf[3] = (uint8_t)(c[3] + 64);
		return 4;
	case 3:
	case 4:
	case 5:
		n = (size_t)1 << logn;
		v = 0;
		xb = 0;
		j = 0;
		for (u = 0; u < n; u += 8) {
			uint32_t lo, hi;

			xb |= encode_ct_8(&lo, &hi, c + u) << j;
			enc24le(buf + v, lo);
			enc32le(buf + v + 3, hi);
			v += 7;
			if (++ j == 8) {
				buf[v ++] = xb;
				xb = 0;
				j = 0;
			}
		}
		if (j != 0) {
			buf[v ++] = xb;
		}
		return v;
	default:
		n = (size_t)1 << logn;
		v = 0;
#if BAT_AVX2
		for (u = 0; u < n; u += 64) {
			buf[v + 56] = encode_ct_8x8(buf + v, c + u);
			v += 57;
		}
#else
		for (u = 0; u < n; u += 64) {
			xb = 0;
			for (j = 0; j < 8; j ++) {
				uint32_t lo, hi;

				xb |= encode_ct_8(&lo, &hi,
					c + u + (j << 3)) << j;
				enc24le(buf + v, lo);
				enc32le(buf + v + 3, hi);
				v += 7;
			}
			buf[v ++] = xb;
		}
#endif
		return v;
	}
}

/* see inner.h */
size_t
bat_decode_ciphertext_257(int8_t *c, unsigned logn,
	const void *in, size_t max_in_len)
{
	const uint8_t *buf;
	uint32_t r, xb;
	size_t in_len, u, n;
#if BAT_AVX2
	__m256i xr;
#endif

	buf = in;
	switch (logn) {
	case 1:
	case 2:
		in_len = (size_t)1 << logn;
		if (max_in_len < in_len) {
			return 0;
		}
		r = 1;
		for (u = 0; u < in_len; u ++) {
			uint8_t tt;

			tt = (uint8_t)(buf[u] - 64);
			c[u] = *(int8_t *)&tt;
			r &= (uint32_t)(buf[u] - 129) >> 31;
		}
		break;
	case 3:
	case 4:
	case 5:
		n = (size_t)1 << logn;
		in_len = 1 + ((size_t)7 << (logn - 3));
		if (max_in_len < in_len) {
			return 0;
		}
		xb = buf[in_len - 1];

		/* Ignored bits in xb must be zero. */
		r = (uint32_t)((xb >> (1u << (logn - 3))) - 1) >> 31;

		/* Decode 57-bit blocks. */
		for (u = 0; u < n; u += 8) {
			size_t v;

			v = u - (u >> 3);
			r &= decode_ct_8(c + u,
				dec24le(buf + v),
				dec32le(buf + v + 3),
				(xb >> ((unsigned)u >> 3)) & 1);
		}
		break;

	default:
		n = (size_t)1 << logn;
		in_len = ((size_t)7 << (logn - 3)) + ((size_t)1 << (logn - 6));
		if (max_in_len < in_len) {
			return 0;
		}
#if BAT_AVX2
		xr = _mm256_set1_epi32(-1);
		for (u = 0; u < n; u += 64) {
			xr = _mm256_and_si256(xr, decode_ct_8x8(c + u, buf));
			buf += 57;
		}
		xr = _mm256_and_si256(xr, _mm256_bsrli_epi128(xr, 4));
		xr = _mm256_and_si256(xr, _mm256_bsrli_epi128(xr, 8));
		r = _mm_cvtsi128_si32(
			_mm_and_si128(
				_mm256_castsi256_si128(xr),
				_mm256_extracti128_si256(xr, 1))) & 1;
#else
		r = 1;
		for (u = 0; u < n; u += 64) {
			size_t v;

			xb = buf[56];
			for (v = 0; v < 64; v += 8) {
				size_t v2;

				v2 = v - (v >> 3);
				r &= decode_ct_8(c + u + v,
					dec24le(buf + v2),
					dec32le(buf + v2 + 3),
					(xb >> ((unsigned)v >> 3)) & 1);
			}
			buf += 57;
		}
#endif
		break;
	}

	return in_len & -(size_t)r;
}
