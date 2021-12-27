#include "inner.h"

/*
 * We use the mod 769 code here.
 */
#define Q   769
#include "modgen.c"

/* see inner.h */
int
bat_make_public_769(uint16_t *h, const int8_t *f, const int8_t *g,
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
bat_encrypt_769(int8_t *c, const uint8_t *sbuf,
	const uint16_t *h, unsigned logn, uint32_t *tmp)
{
	size_t u, n;
	uint16_t *t1, *t2;
	uint32_t e2norm;

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
	 * Coefficients of polynomial s are in {0,1}, extracted from sbuf[].
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
	 * (for q = 769, we have k = 4).
	 *
	 * Rounding is toward +infty, so that error e = k*c - (h*s mod q)
	 * has coefficients in -1..+2 only.
	 *
	 * We also compute in e2norm the sum of the squares of the coefficients
	 * of 2*e', where e' = e - E(e) = e - 1/2*(1+X+X^2+X^3+...+X^(n-1)).
	 */
	e2norm = 0;
	for (u = 0; u < n; u ++) {
		int y, z, ep;

		y = mq_snorm(t1[u]);
		z = ((y + 386) >> 2) - 96;
		c[u] = z;
		ep = (4 * z - y);
		ep = 2 * ep - 1;
		e2norm += (uint32_t)(ep * ep);
	}

	/*
	 * Ciphertext is acceptable if and only if the norm
	 * of (gamma*s',e') is not greater than
	 * 1.08*sqrt((n/2)*gamma^2). Note that e' and s' are
	 * centered on 0, i.e. e' = e - E(e) and s' = s - E(s).
	 * With k = 4, we have gamma = sqrt(5).
	 *
	 * Coefficients of s are in {0,1}, and E(s) = 1/2, therefore
	 * coefficients of s' are in {-1/2,+1/2}. Thus, the sum
	 * of the squares of the coefficients of gamma*s' is always
	 * n*(1/4)*gamma^2 = 5*n/4.
	 *
	 * Coefficients of e are in {-1,0,1,2} and E(e) = 1/2, therefore
	 * coefficients of e' are in {-3/2,-1/2,+1/2,+3/2}. We computed
	 * the sum of the squares of the coefficients of 2*e' in e2norm.
	 *
	 * Thus, the ciphertext is acceptable if and only if:
	 *    5*n/4 + e2norm/4 <= (1.08*sqrt((n/2)*gamma^2))^2
	 * which we rewrite as:
	 *    5*n + e2norm <= 11.664*n
	 * Note that 11.664*n is not an integer for any of the supported
	 * degrees (2 to 1024). For logn in 1..10:
	 *    floor(11.664*2^logn) = floor(11943 / 2^(10-logn))
	 *
	 * We need the comparison to be constant-time, thus we avoid
	 * using a direct comparison operator (which may be translated
	 * by the C compiler into a conditional jump on some architectures).
	 */
	e2norm = (5u << logn) + e2norm - 1u - (11943u >> (10 - logn));
	return e2norm >> 31;
}

/* see inner.h */
void
bat_decrypt_769(uint8_t *sbuf, const int8_t *c,
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
	 * When q = 769, we have Q = 2 and k = 4. This simplifies some
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
	 * possible values for the coefficients of s', which are in
	 * {-1/2,+1/2}. Thus, we can do the computation modulo any prime
	 * p which is not 2, q or q'; in practice, we use the "other q"
	 * (i.e. 257) since we already have the code for computations
	 * modulo that prime.
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
	 * ("Decapsulate"). With q = 769, we have k = 4. Moreover, we
	 * want Q*c, with Q = 2; thus, we compute 8*c here.
	 */
	for (u = 0; u < n; u ++) {
		t1[u] = mq_set(8 * c[u]);
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
		t3[u] = mq_set((int32_t)F[u] * (64513 % 769));
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
		t3[u] = mq_set((int32_t)G[u] * (64513 % 769));
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
	 *  - coefficient must be normalized to -384..+384
	 *  - if the parity is wrong, we must add 769 (if the value is
	 *    negative or zero) or subtract 769 (if the value is strictly
	 *    positive), so that the result is in -768..+769.
	 */
	for (u = 0; u < n; u ++) {
		uint32_t x;

		/*
		 * Get next coefficient of c' into x, normalized in
		 * -384..+384.
		 */
		x = (uint32_t)mq_snorm(t2[u]);

		/*
		 * Adjust x to get the right value modulo 2: if it has
		 * the wrong parity, then we must add or subtract 769.
		 * We add 769 if the value is negative or zero, subtract
		 * 769 if it is negative.
		 */
		x += -(uint32_t)((x ^ cp2) & 1u)
			& -(uint32_t)769
			& (((x - 1) >> 16) & (2 * 769));

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
		 * The Montgomery representation of 1/q' mod q is 602
		 * (with q = 769 and q' = 64513). We add 64596 = 84*769 to
		 * ensure that the value provided to mq_montyred() is in
		 * the proper range (max value will be 602*(768+64596)).
		 */
		x = mq_montyred(602 * (64596 + y0 - y1));

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
		 * If x = 0, set it to q*q' = 49610497.
		 */
		x += 49610497 & -((uint32_t)(x - 1) >> 31);

		/*
		 * Adjust parity to get value modulo 2*q*q': we subtract
		 * q*q' = 49610497 if the value has the wrong parity.
		 */
		x -= 49610497 & -(uint32_t)((x & 1) ^ cs2);

		/*
		 * Value is now in -49610496..+49610497, which is the
		 * correct normalized range. Since we will hand it over
		 * to the module that computes modulo 257, we pre-reduce
		 * it modulo 257. We ensure a positive value by adding
		 * 49610252 = 257 * 193036.
		 */
		t1[u] = m257_tomonty(x + 49610252);
	}

	/*
	 * We now have c' and c'' in t2 and t1, respectively. c' uses
	 * signed integers; c'' is in Montgomery representation modulo 257.
	 * We convert c' the Montgomery representation modulo 257 as
	 * well.
	 */
	for (u = 0; u < n; u ++) {
		t2[u] = m257_tomonty(*(int16_t *)&t2[u] + 257);
	}

	/*
	 * We need to compute:
	 *   Fd = q'*F - f*w
	 *   s' = 1/(q*q'*Q) (Fd*c' - f*c'')
	 *   s = s' + (1/2)*ones
	 *
	 * We use the mod 257 code to obtain q*q'*Q*s'.
	 */
	bat_finish_decapsulate_257(t2, t1, f, F, w, logn, (uint32_t *)t3);

	/*
	 * If the ciphertext is correct and the decapsulation worked well,
	 * then s' has coefficients in {-1/2,1/2} and the coefficients of
	 * s are obtained by adding 1/2. We have the coefficients of
	 * q*q'*Q*s' in t2[], in Montgomery representation modulo 257:
	 *
	 *    s'   s   t2[]
	 *  -1/2   0    12    (-q*q'*Q/2 = 12 mod 257)
	 *  +1/2   1   245    (+q*q'*Q/2 = 245 mod 257)
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
bat_finish_decapsulate_769(uint16_t *cp, uint16_t *cs,
	const int8_t *f, const int8_t *F, const int32_t *w, unsigned logn,
	uint32_t *tmp)
{
	/*
	 * Formulas:
	 *   Fd = q'*F - f*w
	 *   q*q'*Q*s = Fd*c' - f*c''
	 * We use q' = 64513.
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
		 * 87666 = 769 * 114.
		 * This addition ensures that (q' mod q)*F[u] becomes a
		 * positive integer, thus in range for mq_tomonty().
		 */
		t2[u] = mq_tomonty((int32_t)F[u] * (64513 % 769) + 87666);
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
bat_rebuild_G_769(int8_t *G,
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
 * Encoding format: 5 values modulo 769 are encoded over 48 bits.
 * Value x_i is split into its 3 low-order bits (xl_i, value in 0..7)
 * and its 7 high-order bits (xh_i, value in 0..96). First 15 bits
 * are the xl_i; then remaining 33 bits are the xh_i (digits in base-97).
 *
 * Remaining values (0 to 4) are encoded with 10 bits per value.
 */

/*
 * Encode five values x[0]..x[4] into 48 bits (bits 0..31 are returned
 * in *lo, bits 32..47 are the function return value).
 */
static inline uint32_t
encode5(uint32_t *lo, const uint16_t *x)
{
	uint32_t wl, wh, tt, z;
	int i;

	wl = 0;
	for (i = 0; i < 5; i ++) {
		wl |= (uint32_t)(x[i] & 0x07) << (3 * i);
	}

	wh = 0;
	for (i = 4; i > 0; i --) {
		wh = (wh * 97) + (uint32_t)(x[i] >> 3);
	}

	/*
	 * Final iteration may overflow into the extra 'tt' bit.
	 * At that point, wh <= 97^4-1, and 48*wh < 2^32.
	 */
	z = wh + (uint32_t)(x[0] >> 3);
	wh *= 48;
	tt = wh >> 31;
	wh <<= 1;
	wh += z;
	tt |= ((uint32_t)(wh - z) & ~wh) >> 31;

	/*
	 * We have bits 0..14 in wl, bits 15..46 in wh, and bit 47 in tt.
	 */
	*lo = wl | (wh << 15);
	return (wh >> 17) | (tt << 15);
}

#if BAT_AVX2

/*
 * Encode 40 values into 48 bytes (8 times encode5() in parallel).
 */
TARGET_AVX2
static inline void
encode5x8(uint8_t *buf, const uint16_t *x)
{
	__m256i yi0, yi1, yi2, yj0, yj1, yj2, ya0, ya1, ya2, yd0, yd1;
	__m256i yt0, yt1, yt2, yt3, yt4, yt5;
	__m256i yk0, yk1, yk2, yk3, yk4, yk5;
	__m256i ys0, ys1, ys2;
	__m256i y0, y912673w, ym3s, ym16u;

	y0 = _mm256_setzero_si256();
	y912673w = _mm256_set1_epi64x(912673);
	ym3s = _mm256_set1_epi16(0x07);
	ym16u = _mm256_setr_epi16(
		-1, 0, 0, 0, 0, 0, 0, 0,
		-1, 0, 0, 0, 0, 0, 0, 0);

	/*
	 * Load the source values into AVX2 registers.
	 */
	yi0 = _mm256_loadu_si256((const void *)(x +  0));
	yi1 = _mm256_loadu_si256((const void *)(x + 16));
	yi2 = _mm256_castsi128_si256(_mm_loadu_si128((const void *)(x + 32)));

	/*
	 * Extract low 3 bits of each value.
	 */
	ya0 = _mm256_and_si256(yi0, ym3s);
	ya1 = _mm256_and_si256(yi1, ym3s);
	ya2 = _mm256_and_si256(yi2, ym3s);

	/*
	 * Extract high part of each value.
	 */
	yi0 = _mm256_srli_epi16(yi0, 3);
	yi1 = _mm256_srli_epi16(yi1, 3);
	yi2 = _mm256_srli_epi16(yi2, 3);

	/*
	 * Multiply by 97 values with index 1, 3 or 4 modulo 5. This also
	 * clears the unused elements in yi2.
	 */
	yi0 = _mm256_mullo_epi16(yi0, _mm256_setr_epi16(
		1, 97, 1, 97, 97, 1, 97, 1, 97, 97, 1, 97, 1, 97, 97, 1));
	yi1 = _mm256_mullo_epi16(yi1, _mm256_setr_epi16(
		97, 1, 97, 97, 1, 97, 1, 97, 97, 1, 97, 1, 97, 97, 1, 97));
	yi2 = _mm256_mullo_epi16(yi2, _mm256_setr_epi16(
		1, 97, 97, 1, 97, 1, 97, 97, 0, 0, 0, 0, 0, 0, 0, 0));

	/*
	 * Multiply values with index 2 or 3 modulo 5 by 97^2 = 9409.
	 * Results now over 32 bits; we get the low and high 16 bits of
	 * each value in separate registers (yi* for the low parts, yj*
	 * for the high parts).
	 */
	yj0 = _mm256_mulhi_epu16(yi0, _mm256_setr_epi16(
		1, 1, 9409, 9409, 1, 1, 1, 9409,
		9409, 1, 1, 1, 9409, 9409, 1, 1));
	yj1 = _mm256_mulhi_epu16(yi1, _mm256_setr_epi16(
		1, 9409, 9409, 1, 1, 1, 9409, 9409,
		1, 1, 1, 9409, 9409, 1, 1, 1));
	yj2 = _mm256_mulhi_epu16(yi2, _mm256_setr_epi16(
		9409, 9409, 1, 1, 1, 9409, 9409, 1,
		0, 0, 0, 0, 0, 0, 0, 0));
	yi0 = _mm256_mullo_epi16(yi0, _mm256_setr_epi16(
		1, 1, 9409, 9409, 1, 1, 1, 9409,
		9409, 1, 1, 1, 9409, 9409, 1, 1));
	yi1 = _mm256_mullo_epi16(yi1, _mm256_setr_epi16(
		1, 9409, 9409, 1, 1, 1, 9409, 9409,
		1, 1, 1, 9409, 9409, 1, 1, 1));
	yi2 = _mm256_mullo_epi16(yi2, _mm256_setr_epi16(
		9409, 9409, 1, 1, 1, 9409, 9409, 1,
		0, 0, 0, 0, 0, 0, 0, 0));

	/*
	 * Extract the values with index 4 modulo 5, into 64-bit elements
	 * over registers yd0 and yd1, then multiply by 97^3 = 912673
	 * (since these elements have already been multiplied by 97 once).
	 * Each result may use up to 33 bits.
	 */
	yt0 = _mm256_shuffle_epi8(yi0, _mm256_setr_epi8(
		8, 9, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1,
		2, 3, -1, -1, -1, -1, -1, -1,
		12, 13, -1, -1, -1, -1, -1, -1));
	yt1 = _mm256_shuffle_epi8(yi1, _mm256_setr_epi8(
		-1, -1, -1, -1, -1, -1, -1, -1,
		6, 7, -1, -1, -1, -1, -1, -1,
		0, 1, -1, -1, -1, -1, -1, -1,
		10, 11, -1, -1, -1, -1, -1, -1));
	yt2 = _mm256_shuffle_epi8(yi2, _mm256_setr_epi8(
		4, 5, -1, -1, -1, -1, -1, -1,
		14, 15, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1));
	yd0 = _mm256_mul_epu32(y912673w, _mm256_blend_epi32(yt0, yt1, 0x0C));
	yd1 = _mm256_mul_epu32(y912673w, _mm256_blend_epi32(yt1, yt2, 0x0F));

	/*
	 *  yd0:  4  -  -  - 19  -  -  - |  9  -  -  - 14  -  -  -
	 *  yd1: 34  -  -  - 39  -  -  - | 24  -  -  - 29  -  -  -
	 */

	/*
	 * Elements of index 0..3 mod 5 have already been multiplied,
	 * with high parts being in a distinct register. We repack
	 * these values into 32-bit elements.
	 */
	yk0 = _mm256_unpacklo_epi16(yi0, yj0);
	yk1 = _mm256_unpackhi_epi16(yi0, yj0);
	yk2 = _mm256_unpacklo_epi16(yi1, yj1);
	yk3 = _mm256_unpackhi_epi16(yi1, yj1);
	yk4 = _mm256_unpacklo_epi16(yi2, yj2);
	yk5 = _mm256_unpackhi_epi16(yi2, yj2);

	/*
	 *  yk0:  0  1  2  3 |  8  9 10 11
	 *  yk1:  4  5  6  7 | 12 13 14 15
	 *  yk2: 16 17 18 19 | 24 25 26 27
	 *  yk3: 20 21 22 23 | 28 29 30 31
	 *  yk4: 32 33 34 35 |  -  -  -  -
	 *  yk5: 36 37 38 39 |  -  -  -  -
	 */

	/*
	 * Add elements of index 0..3 mod 5 together. Results are
	 * lower than 97^4 = 88529281 < 2^32; we can do so with 32-bit
	 * words. We need to reorder/zero things before adding, so that
	 * each lane contains only values that should indeed be added
	 * together.
	 */

	yt0 = _mm256_blend_epi32(yk0, y0, 0x30);
	yt1 = _mm256_blend_epi32(yk1, y0, 0xC1);
	yt2 = _mm256_blend_epi32(yk2, y0, 0x18);
	yt3 = _mm256_blend_epi32(yk3, y0, 0x30);
	yt4 = _mm256_blend_epi32(yk4, y0, 0x0C);
	yt5 = _mm256_blend_epi32(yk5, y0, 0x08);

	/*
	 *  yt0:  0  1  2  3 |  -  - 10 11
	 *  yt1:  -  5  6  7 | 12 13  -  -
	 *  yt2: 16 17 18  - |  - 25 26 27
	 *  yt3: 20 21 22 23 |  -  - 30 31
	 *  yt4: 32 33  -  - |  -  -  -  -
	 *  yt5: 36 37 38  - |  -  -  -  -
	 */

	yt0 = _mm256_blend_epi32(yt0, yt1, 0x30);
	yt2 = _mm256_blend_epi32(yt2, yk3, 0x10);
	yt5 = _mm256_blend_epi32(yt5, yk4, 0x08);

	/*
	 *  yt0:  0  1  2  3 | 12 13 10 11
	 *  yt1:  -  5  6  7 | 12 13  -  -
	 *  yt2: 16 17 18  - | 28 25 26 27
	 *  yt3: 20 21 22 23 |  -  - 30 31
	 *  yt4: 32 33  -  - |  -  -  -  -
	 *  yt5: 36 37 38 35 |  -  -  -  -
	 */

	yt1 = _mm256_permute4x64_epi64(yt1, 0x4E);
	yk1 = _mm256_permute4x64_epi64(yk1, 0x4E);
	yt4 = _mm256_permute4x64_epi64(yt4, 0x4E);

	/*
	 *  yt0:  0  1  2  3 | 12 13 10 11
	 *  yt1: 12 13  -  - |  -  5  6  7
	 *  yt2: 16 17 18  - | 28 25 26 27
	 *  yt3: 20 21 22 23 |  -  - 30 31
	 *  yt4:  -  -  -  - | 32 33  -  -
	 *  yt5: 36 37 38 35 |  -  -  -  -
	 *
	 *  yk0:  0  1  2  3 |  8  9 10 11
	 *  yk1: 12 13 14 15 |  4  5  6  7
	 *  yk2: 16 17 18 19 | 24 25 26 27
	 *  yk3: 20 21 22 23 | 28 29 30 31
	 *  yk4: 32 33 34 35 |  -  -  -  -
	 *  yk5: 36 37 38 39 |  -  -  -  -
	 */

	yt1 = _mm256_blend_epi32(yt1, yk0, 0x10);
	yt2 = _mm256_blend_epi32(yt2, yk1, 0x08);
	yt3 = _mm256_or_si256(yt3, yt4);
	yt1 = _mm256_blend_epi32(yt1, yt5, 0x0F);

	/*
	 *  yt0:  0  1  2  3 | 12 13 10 11
	 *  yt1: 36 37 38 35 |  8  5  6  7
	 *  yt2: 16 17 18 15 | 28 25 26 27
	 *  yt3: 20 21 22 23 | 32 33 30 31
	 */

	yt0 = _mm256_add_epi32(yt0, _mm256_bsrli_epi128(yt0, 4));
	yt1 = _mm256_add_epi32(yt1, _mm256_bsrli_epi128(yt1, 4));
	yt2 = _mm256_add_epi32(yt2, _mm256_bsrli_epi128(yt2, 4));
	yt3 = _mm256_add_epi32(yt3, _mm256_bsrli_epi128(yt3, 4));
	yt0 = _mm256_add_epi32(yt0, _mm256_bsrli_epi128(yt0, 8));
	yt1 = _mm256_add_epi32(yt1, _mm256_bsrli_epi128(yt1, 8));
	yt2 = _mm256_add_epi32(yt2, _mm256_bsrli_epi128(yt2, 8));
	yt3 = _mm256_add_epi32(yt3, _mm256_bsrli_epi128(yt3, 8));

	/*
	 * We have the sum of elements of index 0..3 mod 5, in the
	 * bottom 32-bit element of lanes of some registers:
	 *
	 *  yt0:  0 | 10
	 *  yt1: 35 |  5
	 *  yt2: 15 | 25
	 *  yt3: 20 | 30
	 *
	 * We now need to promote these values to 64-bit elements,
	 * and dispatch them in the proper emplacements for adding
	 * into yd0 and yd1.
	 */

	yt0 = _mm256_blend_epi32(yt0, y0, 0xEE);
	yt1 = _mm256_blend_epi32(yt1, y0, 0xEE);
	yt2 = _mm256_blend_epi32(yt2, y0, 0xEE);
	yt3 = _mm256_blend_epi32(yt3, y0, 0xEE);

	yt3 = _mm256_permute4x64_epi64(yt3, 0x4E);

	yd0 = _mm256_add_epi64(yd0,
		_mm256_blend_epi32(
			_mm256_or_si256(yt0, _mm256_bslli_epi128(yt2, 8)),
			_mm256_or_si256(yt1, _mm256_bslli_epi128(yt0, 8)),
			0xF0));
	yd1 = _mm256_add_epi64(yd1,
		_mm256_blend_epi32(
			_mm256_or_si256(yt3, _mm256_bslli_epi128(yt1, 8)),
			_mm256_or_si256(yt3, _mm256_bslli_epi128(yt2, 8)),
			0xF0));

	/*
	 * At this point, yd0 and yd1 contain the built 33-bit high parts
	 * of the return values; we need to assemble and adjoin the low
	 * parts. Low bits of the input values are in the ya0, ya1 and ya2
	 * registers (in 16-bit elements).
	 *
	 *  yd0:  0 3 | 1 2
	 *  yd1:  6 7 | 4 5
	 *
	 * To avoid a lot of shuffling and blending, we perform the shifts
	 * with a multiplication.
	 */
	ya0 = _mm256_mullo_epi16(ya0, _mm256_setr_epi16(
		1, 8, 64, 512, 4096, 1, 8, 64,
		512, 4096, 1, 8, 64, 512, 4096, 1));
	ya1 = _mm256_mullo_epi16(ya1, _mm256_setr_epi16(
		8, 64, 512, 4096, 1, 8, 64, 512,
		4096, 1, 8, 64, 512, 4096, 1, 8));
	ya2 = _mm256_mullo_epi16(ya2, _mm256_setr_epi16(
		64, 512, 4096, 1, 8, 64, 512, 4096,
		0, 0, 0, 0, 0, 0, 0, 0));

	/*
	 *  ya0:  0  1  2  3  4  5  6  7 |  8  9 10 11 12 13 14 15
	 *  ya1: 16 17 18 19 20 21 22 23 | 24 25 26 27 28 29 30 31
	 *  ya2: 32 33 34 35 36 37 38 39 |  -  -  -  -  -  -  -  -
	 */

	ys0 = _mm256_setr_epi16(
		-1, -1, -1, -1, -1, 0, 0, 0,
		0, 0, -1, -1, -1, -1, -1, 0);
	yt0 = _mm256_and_si256(ys0, ya0);
	yt1 = _mm256_andnot_si256(ys0, ya0);
	ys1 = _mm256_setr_epi16(
		-1, -1, -1, -1, 0, 0, 0, 0,
		0, -1, -1, -1, -1, -1, 0, 0);
	yt2 = _mm256_and_si256(ys1, ya1);
	yt3 = _mm256_andnot_si256(ys1, ya1);
	ys2 = _mm256_setr_epi16(
		-1, -1, -1, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0);
	yt4 = _mm256_and_si256(ys2, ya2);
	yt5 = _mm256_andnot_si256(ys2, ya2);

	/*
	 *  yt0:  0  1  2  3  4  -  -  - |  -  - 10 11 12 13 14  -
	 *  yt1:  -  -  -  -  -  5  6  7 |  8  9  -  -  -  -  - 15
	 *  yt2: 16 17 18 19  -  -  -  - |  - 25 26 27 28 29  -  -
	 *  yt3:  -  -  -  - 20 21 22 23 | 24  -  -  -  -  - 30 31
	 *  yt4: 32 33 34  -  -  -  -  - |  -  -  -  -  -  -  -  -
	 *  yt5:  -  -  - 35 36 37 38 39 |  -  -  -  -  -  -  -  -
	 */

	yk1 = _mm256_permute4x64_epi64(yt1, 0x4E);
	yk3 = _mm256_permute4x64_epi64(yt3, 0x4E);
	yt5 = _mm256_permute4x64_epi64(yt5, 0x4E);

	yt1 = _mm256_blend_epi32(yt1, yk1, 0xF1);
	yt2 = _mm256_blend_epi32(yt2, yk1, 0x08);
	yt3 = _mm256_blend_epi32(yt3, yk3, 0x01);
	yt4 = _mm256_blend_epi32(yt4, yk3, 0x08);
	yt1 = _mm256_blend_epi32(yt1, yt5, 0xF0);

	/*
	 *  yt0:  0  1  2  3  4  -  -  - |  -  - 10 11 12 13 14  -
	 *  yt1:  8  9  -  -  -  5  6  7 |  -  -  - 35 36 37 38 39
	 *  yt2: 16 17 18 19  -  -  - 15 |  - 25 26 27 28 29  -  -
	 *  yt3: 24  -  -  - 20 21 22 23 | 24  -  -  -  -  - 30 31
	 *  yt4: 32 33 34  -  -  - 30 31 |  -  -  -  -  -  -  -  -
	 */

	yt0 = _mm256_or_si256(yt0, _mm256_srli_epi32(yt0, 16));
	yt1 = _mm256_or_si256(yt1, _mm256_srli_epi32(yt1, 16));
	yt2 = _mm256_or_si256(yt2, _mm256_srli_epi32(yt2, 16));
	yt3 = _mm256_or_si256(yt3, _mm256_srli_epi32(yt3, 16));
	yt4 = _mm256_or_si256(yt4, _mm256_srli_epi32(yt4, 16));
	yt0 = _mm256_or_si256(yt0, _mm256_srli_epi64(yt0, 32));
	yt1 = _mm256_or_si256(yt1, _mm256_srli_epi64(yt1, 32));
	yt2 = _mm256_or_si256(yt2, _mm256_srli_epi64(yt2, 32));
	yt3 = _mm256_or_si256(yt3, _mm256_srli_epi64(yt3, 32));
	yt4 = _mm256_or_si256(yt4, _mm256_srli_epi64(yt4, 32));
	yt0 = _mm256_or_si256(yt0, _mm256_bsrli_epi128(yt0, 8));
	yt1 = _mm256_or_si256(yt1, _mm256_bsrli_epi128(yt1, 8));
	yt2 = _mm256_or_si256(yt2, _mm256_bsrli_epi128(yt2, 8));
	yt3 = _mm256_or_si256(yt3, _mm256_bsrli_epi128(yt3, 8));
	yt4 = _mm256_or_si256(yt4, _mm256_bsrli_epi128(yt4, 8));

	/*
	 * The low parts (15 bits each) are now assembled in the low
	 * 16-bit element of each lane of the yt* registers:
	 *  yt0:  0 | 2
	 *  yt1:  1 | 7
	 *  yt2:  3 | 5
	 *  yt3:  4 | *
	 *  yt4:  6 | -
	 */

	/*
	 * We finally assemble the low and high part of each value,
	 * to get the 48-bit results.
	 *
	 *  yd0:  0 3 | 1 2
	 *  yd1:  6 7 | 4 5
	 */

	yt1 = _mm256_permute4x64_epi64(yt1, 0x4E);
	yt3 = _mm256_permute4x64_epi64(yt3, 0x4E);

	yt0 = _mm256_and_si256(yt0, ym16u);
	yt1 = _mm256_and_si256(yt1, ym16u);
	yt2 = _mm256_and_si256(yt2, ym16u);
	yt3 = _mm256_and_si256(yt3, _mm256_setr_epi16(
		0, 0, 0, 0, 0, 0, 0, 0,
		-1, 0, 0, 0, 0, 0, 0, 0));
	yt4 = _mm256_and_si256(yt4, ym16u);

	/*
	 *  yt0:  0 - | 2 -
	 *  yt1:  7 - | 1 -
	 *  yt2:  3 - | 5 -
	 *  yt3:  - - | 4 -
	 *  yt4:  6 - | - -
	 *
	 *  yd0:  0 3 | 1 2
	 *  yd1:  6 7 | 4 5
	 */

	yd0 = _mm256_slli_epi64(yd0, 15);
	yd1 = _mm256_slli_epi64(yd1, 15);

	yk0 = _mm256_bslli_epi128(yt0, 8);
	yk1 = _mm256_bslli_epi128(yt1, 8);
	yk2 = _mm256_bslli_epi128(yt2, 8);

	yd0 = _mm256_or_si256(yd0,
		_mm256_blend_epi32(
			_mm256_blend_epi32(yt0, yk2, 0x0C),
			_mm256_blend_epi32(yk0, yt1, 0x30), 0xF0));
	yd1 = _mm256_or_si256(yd1,
		_mm256_blend_epi32(
			_mm256_or_si256(yt4, yk1),
			_mm256_or_si256(yt3, yk2), 0xF0));

	/*
	 * We finally have the 48-bit values (in yd0 and yd1). Some
	 * final shuffling will put the bytes in the correct order.
	 *
	 *  yd0:  0 3 | 1 2
	 *  yd1:  6 7 | 4 5
	 */

	yk0 = _mm256_permute4x64_epi64(yd0, 0x4E);
	yk1 = _mm256_permute4x64_epi64(yd1, 0x4E);

	yt0 = _mm256_shuffle_epi8(yd0, _mm256_setr_epi8(
		0, 1, 2, 3, 4, 5, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1,
		12, 13, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1));
	yt1 = _mm256_shuffle_epi8(yk0, _mm256_setr_epi8(
		-1, -1, -1, -1, -1, -1, 0, 1,
		2, 3, 4, 5, 8, 9, 10, 11,
		-1, -1, 8, 9, 10, 11, 12, 13,
		-1, -1, -1, -1, -1, -1, -1, -1));
	yt2 = _mm256_shuffle_epi8(yd1, _mm256_setr_epi8(
		-1, -1, -1, -1, 0, 1, 2, 3,
		4, 5, 8, 9, 10, 11, 12, 13,
		-1, -1, -1, -1, -1, -1, -1, -1,
		0, 1, 2, 3, 4, 5, 8, 9));
	/*
	yt3 = _mm256_shuffle_epi8(yk1, _mm256_setr_epi8(
		10, 11, 12, 13, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1,
		-1, -1, -1, -1, -1, -1, -1, -1));
	*/
	yt3 = _mm256_bsrli_epi128(yk1, 10);

	/*
	 *  yt0:   0  1  2  -  -  -  -  - |  8  -  -  -  -  -  -  -
	 *  yt1:   -  -  -  3  4  5  6  7 |  -  9 10 11  -  -  -  -
	 *  yt2:   -  - 18 19 20 21 22 23 |  -  -  -  - 12 13 14 15
	 *  yt3:  16 17  -  -  -  -  -  - | 22 23  -  -  -  -  -  -
	 */

	yd0 = _mm256_blend_epi32(_mm256_or_si256(yt0, yt1), yt2, 0xC0);
	yd1 = _mm256_or_si256(yt2, yt3);

	_mm256_storeu_si256((void *)(buf +  0), yd0);
	_mm_storeu_si128((void *)(buf + 32), _mm256_castsi256_si128(yd1));
}

#endif

/*
 * Decode a 48-bit word (lo is bits 0..31, hi is bits 32..47)
 * into five values. Returned value is 1 on success, 0 on error.
 */
static uint32_t
decode5(uint16_t *x, uint32_t lo, uint32_t hi)
{
	/*
	 * Useful facts:
	 *    512 = 27 mod 97
	 *    2^19 = 3 mod 97
	 *    for any x in 0..57519, (43241 * x) >> 22 = floor(x / 97)
	 *    97 * 1594008481 = 1 mod 2^32
	 */
	uint32_t a;
	int i;

	/*
	 * First element: we need to compute the quotient and
	 * remainder of the high 33-bit chunk divided by 97.
	 *
	 * We first compute the remainder.
	 */
	a = (lo >> 15) + ((hi & 0x03) << 17)
		+ ((hi >> 2) * 3);             /* 0 <= a <= 573436 */
	a = (a & 0x1FF) + ((a >> 9) * 27);     /* 0 <= a <= 30723 */
	a -= 97 * ((43241 * a) >> 22);
	x[0] = (lo & 0x07) + (a << 3);

	/*
	 * Subtract the remainder from the high 33-bit chunk, and divide by
	 * 97. We can work modulo 2^32 since the division is exact, and the
	 * result fits on 32 bits. From now on, we keep the chunk in 'hi'.
	 */
	hi = (hi << 17) | (lo >> 15);
	hi = (hi - a) * 1594008481u;

	/*
	 * Second element is similar, except that we work on a
	 * 32-bit value now, stored in 'hi' only. Maximum value for 'hi'
	 * is floor((2^33-1) / 97) = 88556026.
	 */
	a = (hi & 0xFFFF) + ((hi >> 16) * 61);  /* 0 <= a <= 147945 */
	a = (a & 0x1FF) + ((a >> 9) * 27);      /* 0 <= a <= 8286 */
	a -= 97 * ((43241 * a) >> 22);
	x[1] = ((lo >> 3) & 0x07) + (a << 3);
	hi = (hi - a) * 1594008481u;

	/*
	 * Third element: maximum value for 'hi' is
	 * floor((2^33-1) / (97^2)) = 912948.
	 */
	a = (hi & 0x0FFF) + ((hi >> 12) * 22);  /* 0 <= a <= 8978 */
	a -= 97 * ((43241 * a) >> 22);
	x[2] = ((lo >> 6) & 0x07) + (a << 3);
	hi = (hi - a) * 1594008481u;

	/*
	 * Fourth element: now, hi <= 9411. We can use a simpler approach.
	 */
	a = (43241 * hi) >> 22;
	x[3] = ((lo >> 9) & 0x07) + ((hi - (97 * a)) << 3);
	hi = a;

	/*
	 * Fifth element: now, hi <= 97. We simply use its value; if it
	 * is off-range, this will be detected in the final test.
	 */
	x[4] = ((lo >> 12) & 0x07) + (hi << 3);

	/*
	 * Decoding is correct if all decoded values are less than 769.
	 */
	a = (uint32_t)-1;
	for (i = 0; i < 5; i ++) {
		a &= (uint32_t)(x[i] - 769);
	}
	return a >> 31;
}

#if BAT_AVX2

/*
 * Decode eight 48-bit words into 40 values. Return value contains eight
 * 32-bit element, with -1 on success and 0 on error.
 */
TARGET_AVX2
static __m256i
decode5x8(uint16_t *x, const uint8_t *src)
{
	__m128i xi0, xi1, xi2;
	__m256i ya0, ya1, yb0, yb1, yc0, yc1, yd0, yd1, yd2, yd3, yd4, yr;
	__m256i y25, y27, y97, y769, y43241, y1594008481;
	__m256i ym3, ym3w, ym9, ym17, ym19w;

	y25 = _mm256_set1_epi32(25);
	y27 = _mm256_set1_epi32(27);
	y97 = _mm256_set1_epi32(97);
	y769 = _mm256_set1_epi32(769);
	y43241 = _mm256_set1_epi32(43241);
	y1594008481 = _mm256_set1_epi32(1594008481);
	ym3 = _mm256_set1_epi32(0x07);
	ym3w = _mm256_set1_epi64x(0x07);
	ym9 = _mm256_set1_epi32(0x1FF);
	ym17 = _mm256_set1_epi32(0x1FFFF);
	ym19w = _mm256_set1_epi64x(0x7FFFF);

	/*
	 * Read the next 48 bytes.
	 */
	xi0 = _mm_loadu_si128((const void *)(src +  0));
	xi1 = _mm_loadu_si128((const void *)(src + 16));
	xi2 = _mm_loadu_si128((const void *)(src + 32));

	/*
	 * Useful facts:
	 *    512 = 27 mod 97
	 *    2^19 = 3 mod 97
	 *    for any x in 0..57519, (43241 * x) >> 22 = floor(x / 97)
	 *    97 * 1594008481 = 1 mod 2^32
	 */

	/*
	 * Split the eight input values over eight 64-bit elements
	 * (into ya0 and ya1). We also extract the 33-bit high part of
	 * each 48-bit input chunk (into yb0 and yb1).
	 */
	ya0 = _mm256_setr_m128i(
		_mm_shuffle_epi8(xi0, _mm_setr_epi8(
			0, 1, 2, 3, 4, 5, -1, -1,
			6, 7, 8, 9, 10, 11, -1, -1)),
		_mm_or_si128(
			_mm_bsrli_si128(xi0, 12),
			_mm_shuffle_epi8(xi1, _mm_setr_epi8(
				-1, -1, -1, -1, 0, 1, -1, -1,
				2, 3, 4, 5, 6, 7, -1, -1))));
	ya1 = _mm256_setr_m128i(
		_mm_or_si128(
			_mm_shuffle_epi8(xi1, _mm_setr_epi8(
				8, 9, 10, 11, 12, 13, -1, -1,
				14, 15, -1, -1, -1, -1, -1, -1)),
			_mm_shuffle_epi8(xi2, _mm_setr_epi8(
				-1, -1, -1, -1, -1, -1, -1, -1,
				-1, -1, 0, 1, 2, 3, -1, -1))),
		_mm_shuffle_epi8(xi2, _mm_setr_epi8(
			4, 5, 6, 7, 8, 9, -1, -1,
			10, 11, 12, 13, 14, 15, -1, -1)));
	yb0 = _mm256_srli_epi64(ya0, 15);
	yb1 = _mm256_srli_epi64(ya1, 15);

	/*
	 * Compute the remainder (mod 97) of each 33-bit high part.
	 *
	 *   x <- (x mod 2^19) + 3*floor(x / 2^19)    (x <= 573436)
	 *   x <- (x mod 2^9) + 27*floor(x / 2^9)     (x <= 30723)
	 *   x <- x - 97 * floor((43241 * x) / 2^22)
	 */
	yc0 = _mm256_srli_epi64(yb0, 19);
	yc1 = _mm256_srli_epi64(yb1, 19);
	yc0 = _mm256_add_epi64(
		_mm256_add_epi64(_mm256_and_si256(yb0, ym19w), yc0),
		_mm256_add_epi64(yc0, yc0));
	yc1 = _mm256_add_epi64(
		_mm256_add_epi64(_mm256_and_si256(yb1, ym19w), yc1),
		_mm256_add_epi64(yc1, yc1));
	yc0 = _mm256_add_epi64(
		_mm256_and_si256(yc0, ym9),
		_mm256_mullo_epi16(_mm256_srli_epi64(yc0, 9), y27));
	yc1 = _mm256_add_epi64(
		_mm256_and_si256(yc1, ym9),
		_mm256_mullo_epi16(_mm256_srli_epi64(yc1, 9), y27));
	yc0 = _mm256_sub_epi64(yc0,
		_mm256_mullo_epi16(y97,
			_mm256_srli_epi64(_mm256_mulhi_epu16(yc0, y43241), 6)));
	yc1 = _mm256_sub_epi64(yc1,
		_mm256_mullo_epi16(y97,
			_mm256_srli_epi64(_mm256_mulhi_epu16(yc1, y43241), 6)));

	/*
	 * Get first element of each group of 5, repacked into 32-bit
	 * elements.
	 *
	 *   yd0:  e0 e4 e1 e5 | e2 e6 e3 e7
	 */
	yd0 = _mm256_or_si256(
		_mm256_or_si256(
			_mm256_and_si256(ya0, ym3w),
			_mm256_slli_epi32(yc0, 3)),
		_mm256_bslli_epi128(
			_mm256_or_si256(
				_mm256_and_si256(ya1, ym3w),
				_mm256_slli_epi32(yc1, 3)), 4));

	/*
	 * Subtract the remainder from the high 33-bit value and
	 * divide by 97, which we can do with a multiplication since
	 * the division is exact. Moreover, the result will fit on
	 * 32 bits so we can work modulo 2^32, i.e. repack into a
	 * single AVX2 register (yb0). We also repack the low 15-bit
	 * words (ya0).
	 */
	yb0 = _mm256_mul_epu32(_mm256_sub_epi64(yb0, yc0), y1594008481);
	yb1 = _mm256_mul_epu32(_mm256_sub_epi64(yb1, yc1), y1594008481);
	ya0 = _mm256_blend_epi32(ya0, _mm256_bslli_epi128(ya1, 4), 0xAA);
	yb0 = _mm256_blend_epi32(yb0, _mm256_bslli_epi128(yb1, 4), 0xAA);

	/*
	 * We now use 32-bit elements, in the same order as yd0. Maximum
	 * value of each element is floor((2^33-1) / 97) = 88556026.
	 *
	 *   x <- (x mod 2^17) + 25*floor(x / 2^17)    (x <= 147921)
	 *   x <- (x mod 2^9) + 27*floor(x / 2^9)      (x <= 8260)
	 *   x <- x - 97 * floor((43241 * x) / 2^22)
	 *
	 * Note that in '25*floor(x / 2^17)', the second operand is at
	 * most 675, so the multiplication result fits on 16 bits.
	 * Similarly, in '27*floor(x / 2^9)', the second operand is at
	 * most 288, and the multiplication can again be done on 16-bit
	 * values.
	 *
	 * We produce the second element of each group in yd1, while we
	 * also subtract the remainder from the current values and divide
	 * by 97.
	 */
	yc0 = _mm256_add_epi32(
		_mm256_and_si256(yb0, ym17),
		_mm256_mullo_epi16(_mm256_srli_epi32(yb0, 17), y25));
	yc0 = _mm256_add_epi32(
		_mm256_and_si256(yc0, ym9),
		_mm256_mullo_epi16(_mm256_srli_epi32(yc0, 9), y27));
	yc0 = _mm256_sub_epi32(yc0,
		_mm256_mullo_epi16(y97,
			_mm256_srli_epi32(_mm256_mulhi_epu16(yc0, y43241), 6)));
	yb0 = _mm256_mullo_epi32(_mm256_sub_epi32(yb0, yc0), y1594008481);
	yd1 = _mm256_or_si256(
		_mm256_and_si256(_mm256_srli_epi32(ya0, 3), ym3),
		_mm256_slli_epi32(yc0, 3));

	/*
	 * Third element is produced with the same method. This time, the
	 * current value is at most floor((2^33-1) / 97^2) = 912948.
	 */
	yc0 = _mm256_add_epi32(
		_mm256_and_si256(yb0, ym17),
		_mm256_mullo_epi16(_mm256_srli_epi32(yb0, 17), y25));
	yc0 = _mm256_add_epi32(
		_mm256_and_si256(yc0, ym9),
		_mm256_mullo_epi16(_mm256_srli_epi32(yc0, 9), y27));
	yc0 = _mm256_sub_epi32(yc0,
		_mm256_mullo_epi16(y97,
			_mm256_srli_epi32(_mm256_mulhi_epu16(yc0, y43241), 6)));
	yb0 = _mm256_mullo_epi32(_mm256_sub_epi32(yb0, yc0), y1594008481);
	yd2 = _mm256_or_si256(
		_mm256_and_si256(_mm256_srli_epi32(ya0, 6), ym3),
		_mm256_slli_epi32(yc0, 3));

	/*
	 * Input value is now at most floor((2^33 - 1) / 97^3) = 9411,
	 * so we can use a simpler approach for the 4th and 5th elements.
	 */
	yc0 = _mm256_srli_epi32(_mm256_mulhi_epu16(yb0, y43241), 6);
	yc1 = _mm256_sub_epi32(yb0, _mm256_mullo_epi16(y97, yc0));
	yd4 = _mm256_or_si256(
		_mm256_and_si256(_mm256_srli_epi32(ya0, 12), ym3),
		_mm256_slli_epi32(yc0, 3));
	yd3 = _mm256_or_si256(
		_mm256_and_si256(_mm256_srli_epi32(ya0, 9), ym3),
		_mm256_slli_epi32(yc1, 3));

	/*
	 * Compute status result: all output elements should be less
	 * than 769.
	 */
	yr = _mm256_and_si256(
		_mm256_cmpgt_epi32(y769, yd0),
		_mm256_and_si256(
			_mm256_and_si256(
				_mm256_cmpgt_epi32(y769, yd1),
				_mm256_cmpgt_epi32(y769, yd2)),
			_mm256_and_si256(
				_mm256_cmpgt_epi32(y769, yd3),
				_mm256_cmpgt_epi32(y769, yd4))));

	/*
	 * We now have all values to return, in yd0..yd4 (32 bits per
	 * element):
	 *
	 *   yd0:  0 20  5 25 | 10 30 15 35
	 *   yd1:  1 21  6 26 | 11 31 16 36
	 *   yd2:  2 22  7 27 | 12 32 17 37
	 *   yd3:  3 23  8 28 | 13 33 18 38
	 *   yd4:  4 24  9 29 | 14 34 19 39
	 *
	 * We need to repack them in the proper order, with 16 bits per
	 * element. We first reassemble pairs, which reduces the number
	 * of registers to handle.
	 */
	yc0 = _mm256_or_si256(yd0, _mm256_slli_epi32(yd1, 16));
	yc1 = _mm256_or_si256(yd2, _mm256_slli_epi32(yd3, 16));

	/*
	 * yc0:   0  1 20 21  5  6 25 26 | 10 11 30 31 15 16 35 36
	 * yc1:   2  3 22 23  7  8 27 28 | 12 13 32 33 17 18 37 38
	 * yd4:   4  - 24  -  9  - 29  - | 14  - 34  - 19  - 39  -  
	 */

	yd0 = _mm256_or_si256(
		_mm256_or_si256(
			_mm256_shuffle_epi8(yc0, _mm256_setr_epi8(
				0, 1, 2, 3, -1, -1, -1, -1,
				-1, -1, 8, 9, 10, 11, -1, -1,
				-1, -1, -1, -1, 0, 1, 2, 3,
				-1, -1, -1, -1, -1, -1, 8, 9)),
			_mm256_shuffle_epi8(yc1, _mm256_setr_epi8(
				-1, -1, -1, -1, 0, 1, 2, 3,
				-1, -1, -1, -1, -1, -1, 8, 9,
				-1, -1, -1, -1, -1, -1, -1, -1,
				0, 1, 2, 3, -1, -1, -1, -1))),
		_mm256_shuffle_epi8(yd4, _mm256_setr_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1,
			0, 1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, 0, 1, -1, -1)));
	yd1 = _mm256_or_si256(
		_mm256_or_si256(
			_mm256_shuffle_epi8(yc0, _mm256_setr_epi8(
				-1, -1, -1, -1, -1, -1, -1, -1,
				4, 5, 6, 7, -1, -1, -1, -1,
				10, 11, -1, -1, -1, -1, -1, -1,
				-1, -1, -1, -1, 4, 5, 6, 7)),
			_mm256_shuffle_epi8(yc1, _mm256_setr_epi8(
				10, 11, -1, -1, -1, -1, -1, -1,
				-1, -1, -1, -1, 4, 5, 6, 7,
				-1, -1, 8, 9, 10, 11, -1, -1,
				-1, -1, -1, -1, -1, -1, -1, -1))),
		_mm256_shuffle_epi8(yd4, _mm256_setr_epi8(
			-1, -1, 8, 9, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, 8, 9,
			-1, -1, -1, -1, -1, -1, -1, -1)));
	yd2 = _mm256_or_si256(
		_mm256_or_si256(
			_mm256_shuffle_epi8(yc0, _mm256_setr_epi8(
				-1, -1, 12, 13, 14, 15, -1, -1,
				-1, -1, -1, -1, -1, -1, -1, -1,
				-1, -1, -1, -1, -1, -1, 12, 13,
				14, 15, -1, -1, -1, -1, -1, -1)),
			_mm256_shuffle_epi8(yc1, _mm256_setr_epi8(
				-1, -1, -1, -1, -1, -1, 12, 13,
				14, 15, -1, -1, -1, -1, -1, -1,
				4, 5, 6, 7, -1, -1, -1, -1,
				-1, -1, 12, 13, 14, 15, -1, -1))),
		_mm256_shuffle_epi8(yd4, _mm256_setr_epi8(
			4, 5, -1, -1, -1, -1, -1, -1,
			-1, -1, 12, 13, -1, -1, -1, -1,
			-1, -1, -1, -1, 4, 5, -1, -1,
			-1, -1, -1, -1, -1, -1, 12, 13)));

	/*
	 * yd0:   0  1  2  3  4  5  6  7 |  -  - 10 11 12 13 14 15
	 * yd1:   8  9  -  - 20 21 22 23 | 16 17 18 19  -  - 30 31
	 * yd2:  24 25 26 27 28 29  -  - | 32 33 34 35 36 37 38 39
	 */

	yc0 = _mm256_permute4x64_epi64(yd1, 0x4E);
	yd2 = _mm256_permute4x64_epi64(yd2, 0x4E);

	/*
	 * yd0:   0  1  2  3  4  5  6  7 |  -  - 10 11 12 13 14 15
	 * yd1:   8  9  -  - 20 21 22 23 | 16 17 18 19  -  - 30 31
	 * yd2:  32 33 34 35 36 37 38 39 | 24 25 26 27 28 29  -  -
	 * yc0:  16 17 18 19  -  - 30 31 |  8  9  -  - 20 21 22 23
	 */

	yd0 = _mm256_blend_epi32(yd0, yc0, 0x10);
	yd1 = _mm256_blend_epi32(
		_mm256_blend_epi32(yd1, yc0, 0x03),
		yd2, 0x70);

	_mm256_storeu_si256((void *)(x +  0), yd0);
	_mm256_storeu_si256((void *)(x + 16), yd1);
	_mm_storeu_si128((void *)(x + 32), _mm256_castsi256_si128(yd2));
	return yr;
}

#endif

/* see inner.h */
size_t
bat_encode_769(void *out, size_t max_out_len,
	const uint16_t *x, unsigned logn)
{
	size_t out_len, u, v, n;
	uint8_t *buf;
	uint32_t acc;
	int acc_len;

	n = (size_t)1 << logn;
	out_len = 6 * (n / 5) + (n % 5) + 1;
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
	default:
		v = 0;
		u = 0;
#if BAT_AVX2
		for (; (u + 40) <= n; u += 40) {
			encode5x8(buf + v, x + u);
			v += 48;
		}
#endif
		for (; (u + 5) <= n; u += 5) {
			uint32_t lo, hi;

			hi = encode5(&lo, x + u);
			enc32le(buf + v, lo);
			enc16le(buf + v + 4, hi);
			v += 6;
		}
		acc = 0;
		acc_len = 0;
		while (u < n) {
			acc |= (uint32_t)x[u ++] << acc_len;
			buf[v ++] = (uint8_t)acc;
			acc >>= 8;
			acc_len += 2;
		}
		if (acc_len != 0) {
			buf[v ++] = (uint8_t)acc;
		}
		return v;
	}
}

/* see inner.h */
size_t
bat_decode_769(uint16_t *x, unsigned logn,
	const void *in, size_t max_in_len)
{
	const uint8_t *buf;
	uint32_t r;
	size_t in_len, u, n;
#if BAT_AVX2
	__m256i yr;
#endif

	buf = in;
	n = (size_t)1 << logn;
	in_len = 6 * (n / 5) + (n % 5) + 1;
	if (max_in_len < in_len) {
		return 0;
	}
	u = 0;
#if BAT_AVX2
	yr = _mm256_set1_epi32(-1);
	for (; (u + 40) <= n; u += 40) {
		yr = _mm256_and_si256(yr, decode5x8(x + u, buf));
		buf += 48;
	}
	yr = _mm256_and_si256(yr, _mm256_bsrli_epi128(yr, 4));
	yr = _mm256_and_si256(yr, _mm256_bsrli_epi128(yr, 8));
	r = _mm_cvtsi128_si32(
		_mm_and_si128(
			_mm256_castsi256_si128(yr),
			_mm256_extracti128_si256(yr, 1))) & 1;
#else
	r = 1;
#endif
	for (; (u + 5) <= n; u += 5) {
		r &= decode5(x + u, dec32le(buf), dec16le(buf + 4));
		buf += 6;
	}
	if (u < n) {
		uint32_t acc;
		int acc_len;

		acc = *buf ++;
		acc_len = 8;
		while (u < n) {
			uint32_t z;

			acc |= (uint32_t)(*buf ++) << acc_len;
			z = acc & 0x3FF;
			acc >>= 10;
			acc_len -= 2;
			x[u ++] = z;
			r &= (uint32_t)(z - 769) >> 31;
		}
		r &= (uint32_t)(acc - 1) >> 31;
	}

	return in_len & -(size_t)r;
}

/*
 * Ciphertext encoding format: 5 values are encoded over 38 bits.
 * Each value is in -96..+96; we first add 96 to get a value in 0..192.
 * Values are then encoded in base 193; 193^5 is slightly below 2^38.
 *
 * Remaining values (0 to 4) are encoded with 8 bits per value.
 */

/*
 * Encode five values x[0]..x[4] into 38 bits (bits 0..31 are returned
 * in *lo, bits 32..37 are the function return value).
 */
static inline uint32_t
encode_ct_5(uint32_t *lo, const int8_t *c)
{
	uint32_t wl, wh, x, y;

	/*
	 * First four values fit in 30 bits.
	 */
	wl = (uint32_t)(c[0] + 96)
		+ 193 * (uint32_t)(c[1] + 96)
		+ (193 * 193) * (uint32_t)(c[2] + 96)
		+ (193 * 193 * 193) * (uint32_t)(c[3] + 96);

	/*
	 * We want to add 193^4 * x (with x = c[4] + 96, in 0..192).
	 * Multiplication result may be above 2^32, so we need to split
	 * things in order to remain on pure 32-bit code. Note that:
	 *
	 *   193^4 = 0x52B36301 = (0x52B363 * 256) + 1
	 *
	 * Therefore, for 0 <= x <= 192:
	 *
	 *   x * (193^4) = (0x52B363 * x) * 256 + x
	 *
	 * Note that 0x52B363 * x will be less than 2^30, and the final
	 * addition of x cannot propagate any carry into bits 8+.
	 */
	x = (uint32_t)(c[4] + 96);
	y = 0x52B363 * x;
	wh = y >> 24;
	y = (y << 8) + x;

	/*
	 * Addition of y with wl may trigger a carry. Note that wl, at
	 * this point, is at most 193^4 - 1, which is lower than 2^30;
	 * thus, if a carry occurs, then the result will be less than
	 * 2^30 as well. We know a carry occured if and only if bit 31
	 * of y is zero AND bit 31 of y - wl is 1.
	 */
	y += wl;
	wh += (uint32_t)((y - wl) & ~y) >> 31;

	*lo = y;
	return wh;
}

#if BAT_AVX2

/*
 * Encode 40 values x[0]..x[39] into 8*38 = 304 bits = 38 bytes.
 */
static inline void
encode_ct_5x8(uint8_t *dst, const int8_t *c)
{
	__m256i xi0, xi1, xv0, xv1, xv2, xv3, xv4, xv5;
	__m256i x96, x193, x37249, x7189057;

	x96 = _mm256_set1_epi8(96);
	x193 = _mm256_set1_epi32(193);
	x37249 = _mm256_set1_epi32(37249);
	x7189057 = _mm256_set1_epi32(7189057);

	/*
	 * Load the 40 values into two registers and normalize to 0..192
	 * the values.
	 */
	xi0 = _mm256_loadu_si256((const void *)c);
	xi1 = _mm256_castsi128_si256(_mm_loadu_si128((const void *)(c + 24)));
	xi0 = _mm256_add_epi8(xi0, x96);
	xi1 = _mm256_add_epi8(xi1, x96);
	xi1 = _mm256_permute2x128_si256(xi0, xi1, 0x21);

	/*
	 * xi0:  0 .. 15 | 16 .. 31
	 * xi1: 16 .. 31 | 24 .. 39
	 */

	/*
	 * Split the 40 values over 5 AVX2 registers, 32 bits per value.
	 */
	xv0 = _mm256_or_si256(
		_mm256_shuffle_epi8(xi0, _mm256_setr_epi8(
			0, -1, -1, -1, 5, -1, -1, -1,
			10, -1, -1, -1, 15, -1, -1, -1,
			4, -1, -1, -1, 9, -1, -1, -1,
			14, -1, -1, -1, -1, -1, -1, -1)),
		_mm256_shuffle_epi8(xi1, _mm256_setr_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, 11, -1, -1, -1)));
	xv1 = _mm256_or_si256(
		_mm256_shuffle_epi8(xi0, _mm256_setr_epi8(
			1, -1, -1, -1, 6, -1, -1, -1,
			11, -1, -1, -1, -1, -1, -1, -1,
			5, -1, -1, -1, 10, -1, -1, -1,
			15, -1, -1, -1, -1, -1, -1, -1)),
		_mm256_shuffle_epi8(xi1, _mm256_setr_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, 0, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, 12, -1, -1, -1)));
	xv2 = _mm256_or_si256(
		_mm256_shuffle_epi8(xi0, _mm256_setr_epi8(
			2, -1, -1, -1, 7, -1, -1, -1,
			12, -1, -1, -1, -1, -1, -1, -1,
			6, -1, -1, -1, 11, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1)),
		_mm256_shuffle_epi8(xi1, _mm256_setr_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, 1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1,
			8, -1, -1, -1, 13, -1, -1, -1)));
	xv3 = _mm256_or_si256(
		_mm256_shuffle_epi8(xi0, _mm256_setr_epi8(
			3, -1, -1, -1, 8, -1, -1, -1,
			13, -1, -1, -1, -1, -1, -1, -1,
			7, -1, -1, -1, 12, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1)),
		_mm256_shuffle_epi8(xi1, _mm256_setr_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, 2, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1,
			9, -1, -1, -1, 14, -1, -1, -1)));
	xv4 = _mm256_or_si256(
		_mm256_shuffle_epi8(xi0, _mm256_setr_epi8(
			4, -1, -1, -1, 9, -1, -1, -1,
			14, -1, -1, -1, -1, -1, -1, -1,
			8, -1, -1, -1, 13, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1)),
		_mm256_shuffle_epi8(xi1, _mm256_setr_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, 3, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1,
			10, -1, -1, -1, 15, -1, -1, -1)));

	/*
	 * Do multiplications:
	 *   xv0 is unchanged
	 *   xv1 <- 193*xv1
	 *   xv2 <- 193*193*xv2
	 *   xv3 <- 193*193*193*xv3
	 *   xv4 <- 193*193*193*193*xv4
	 * Result of 193*xv1 fits on 16 bits, so we can use mullo_epi16()
	 * which has lower latency.
	 * For xv4, result does not fit into 32 bits, so we must do it
	 * differently. We first multiply it with 193^3 instead of 193^4.
	 */
	xv1 = _mm256_mullo_epi16(xv1, x193);
	xv2 = _mm256_mullo_epi32(xv2, x37249);
	xv3 = _mm256_mullo_epi32(xv3, x7189057);
	xv4 = _mm256_mullo_epi32(xv4, x7189057);

	/*
	 * 193 = 1 + 3*(2^6) so the result should be:
	 *   xv0 + xv1 + xv2 + xv3 + xv4 + 3*(xv4 << 6)
	 * Note that max(xv0 + xv1 + xv2 + xv3 + xv4) =
	 * 192 + 192*193 + 192*193^2 + 2*192*193^3 = 2767786944, which is
	 * lower than 2^32. The low 6 bits of this value are moreover
	 * correct since 3*(xv4 << 6) is a multiple of 64. We can thus
	 * get the low 6 bits separately, and shift to make room for the
	 * final additions.
	 */
	xv0 = _mm256_add_epi32(
		xv4,
		_mm256_add_epi32(
			_mm256_add_epi32(xv0, xv1),
			_mm256_add_epi32(xv2, xv3)));
	xv4 = _mm256_add_epi32(xv4, _mm256_add_epi32(xv4, xv4));
	xv1 = _mm256_add_epi32(
		_mm256_srli_epi32(xv0, 6), xv4);
	xv0 = _mm256_add_epi32(
		_mm256_slli_epi32(xv4, 6), xv0);

	/*
	 * Get some shifted versions:
	 *
	 *   xv0: bits 0..31
	 *   xv1: bits 6..37
	 *   xv2: bits 16..37 | bits 0..9 of next value
	 *   xv3: bits 10..37 | bits 0..3 of next value
	 *   xv4: bits 4..35
	 *   xv5: bits 36..37 of previous value | bits 0..29
	 */
	xv2 = _mm256_or_si256(
		_mm256_srli_epi32(xv1, 10),
		_mm256_bsrli_epi128(_mm256_slli_epi32(xv0, 22), 4));
	xv3 = _mm256_or_si256(
		_mm256_srli_epi32(xv1, 4),
		_mm256_bsrli_epi128(_mm256_slli_epi32(xv0, 28), 4));
	xv4 = _mm256_or_si256(
		_mm256_srli_epi32(xv0, 4),
		_mm256_slli_epi32(xv1, 2));
	xv5 = _mm256_or_si256(
		_mm256_slli_epi32(xv0, 2),
		_mm256_bslli_epi128(_mm256_srli_epi32(xv1, 30), 4));

	/*
	 * We have our eight 38-bit values, currenly over xv0 and xv1;
	 * xv0 contains bits 0..31 and xv1 contains bits 6..37. We
	 * assemble in xi0 and xi1 the two 19-byte halves of output
	 * (low lanes contain the first chunk, high lanes the second chunk).
	 *
	 * Output byte provenance (low lane):
	 *    0: xv0[0]
	 *    1: xv0[1]
	 *    2: xv0[2]
	 *    3: xv0[3]
	 *    4: xv2[2]
	 *    5: xv2[3]
	 *    6: xv3[4]
	 *    7: xv3[5]
	 *    8: xv3[6]
	 *    9: xv3[7]
	 *   10: xv4[8]
	 *   11: xv4[9]
	 *   12: xv4[10]
	 *   13: xv4[11]
	 *   14: xv5[12]
	 *   15: xv5[13]
	 *   16: xv1[13]
	 *   17: xv1[14]
	 *   18: xv1[15]
	 */
	xi0 = _mm256_blend_epi32(
		xv0,
		_mm256_bslli_epi128(
			_mm256_blend_epi32(
				_mm256_blend_epi32(xv2, xv3, 0xEE),
				_mm256_blend_epi32(xv4, xv5, 0x88), 0xCC), 2),
		0xEE);
	xi1 = _mm256_bsrli_epi128(xv1, 13);

	/*
	 * Assemble all the output bytes:
	 *   16 bytes of low lane of xi0
	 *   3 bytes of low lane of xi1
	 *   16 bytes of high lane of xi0
	 *   3 bytes of high lane of xi1
	 */

	/* xv2: 19..34 | 16..18 0-0 */
	xv2 = _mm256_permute2x128_si256(xi0, xi1, 0x21);

	/* xv3: 35..37 0-0 | 35..37 0-0 */
	xv3 = _mm256_permute4x64_epi64(xi1, 0xEE);

	/* xv0:  0..15 | 16..31 */
	xv0 = _mm256_blend_epi32(
		xi0,
		_mm256_or_si256(_mm256_bslli_epi128(xi0, 3), xv2),
		0xF0);

	/* xv1: 32..37 undef | undef */
	xv1 = _mm256_or_si256(
		_mm256_bslli_epi128(xv3, 3),
		_mm256_bsrli_epi128(xv2, 13));

	/*
	 * Write results.
	 */
	_mm256_storeu_si256((void *)dst, xv0);
	*(uint32_t *)(dst + 32) =
		_mm_cvtsi128_si32(_mm256_castsi256_si128(xv1));
	*(uint16_t *)(dst + 36) =
		_mm_cvtsi128_si32(_mm256_castsi256_si128(
			_mm256_bsrli_epi128(xv1, 4)));
}

#endif

/*
 * Decode a 38-bit word (lo is bits 0..31, hi is bits 32..37)
 * into five values. Returned value is 1 on success, 0 on error.
 */
static uint32_t
decode_ct_5(int8_t *c, uint32_t lo, uint32_t hi)
{
	/*
	 * Useful facts:
	 *
	 *   2^8 = 63 mod 193
	 *   2^16 = 109 mod 193
	 *   2^24 = 112 mod 193
	 *   2^32 = 108 mod 193
	 *   (x * 43465) >> 23 == floor(x / 193)  for all x in 0..61372.
	 *
	 * General processing is the following:
	 *
	 *  - Input is an integer z.
	 *  - We find a = z mod 193, which yields the next ciphertext value.
	 *  - z <- (z - a) / 193. This is an exact division, and the result
	 *    fits on 32 bits, so we can compute it by multiplying z
	 *    with 2425655105, which is the inverse of 193 modulo 2^32.
	 *    We can even do that by working with z mod 2^32.
	 */
	uint32_t a;

	/*
	 * First iteration: lo + 2^32*hi <= 2^38-1
	 */
	a = (lo & 0xFF)
		+ 63 * ((lo >> 8) & 0xFF)
		+ 109 * ((lo >> 16) & 0xFF)
		+ 112 * ((lo >> 24) & 0xFF)
		+ 108 * hi;                    /* a = 0..79479 */
	a = (a & 0x3FF) + 59 * (a >> 10);      /* a = 0..5565 */
	a -= 193 * ((a * 43465) >> 23);
	c[0] = (int)a - 96;
	lo = (lo - a) * 2425655105u;

	/*
	 * Second iteration: lo <= floor((2^38-1) / 193) = 1424237859
	 */
	a = (lo & 0xFF)
		+ 63 * ((lo >> 8) & 0xFF)
		+ 109 * ((lo >> 16) & 0xFF)
		+ 112 * ((lo >> 24) & 0xFF);   /* a = 0..53522 */
	a -= 193 * ((a * 43465) >> 23);
	c[1] = (int)a - 96;
	lo = (lo - a) * 2425655105u;

	/*
	 * Third iteration: lo <= floor((2^38-1) / 193^2) = 7379470
	 */
	a = (lo & 0xFF)
		+ 63 * ((lo >> 8) & 0xFF)
		+ 109 * ((lo >> 16) & 0xFF);   /* a = 0..28527 */
	a -= 193 * ((a * 43465) >> 23);
	c[2] = (int)a - 96;
	lo = (lo - a) * 2425655105u;

	/*
	 * Fourth iteration: lo <= floor((2^38-1) / 193^3) = 38235
	 */
	a = (lo * 43465) >> 23;
	c[3] = (int)(lo - 193 * a) - 96;
	lo = a;

	/*
	 * Fifth iteration: lo <= 198.
	 * Decoding is correct if and only if lo <= 192 at this point.
	 */
	c[4] = (int)lo - 96;
	return (uint32_t)(lo - 193) >> 31;
}

#if BAT_AVX2

/*
 * Decode eight 38-bit words in parallel (38 input bytes), yielding
 * 40 ciphertext elements. Returned value is a set of eight 32-bit
 * elements (-1 on success, 0 on failure).
 */
TARGET_AVX2
static __m256i
decode_ct_5x8(int8_t *c, const uint8_t *src)
{
	__m256i yi0, yi1, ya0, ya1, yd0, yd1, yd2, yd3, yd4, yr;
	__m256i y59, y63, y96n, y108, y109, y112, y193, y43465, y2425655105;
	__m256i ym8, ym8w, ym10, ym38w;

	y59 = _mm256_set1_epi32(59);
	y63 = _mm256_set1_epi32(63);
	y96n = _mm256_set1_epi8(96);
	y108 = _mm256_set1_epi32(108);
	y109 = _mm256_set1_epi32(109);
	y112 = _mm256_set1_epi32(112);
	y193 = _mm256_set1_epi32(193);
	y43465 = _mm256_set1_epi32(43465);
	y2425655105 = _mm256_set1_epi32(-1869312191);
	ym8 = _mm256_set1_epi32(0xFF);
	ym8w = _mm256_set1_epi64x(0xFF);
	ym10 = _mm256_set1_epi32(0x3FF);
	ym38w = _mm256_set1_epi64x(0x3FFFFFFFFF);

	/*
	 * Load all bytes and split into 38-bit words, each in a
	 * 64-bit element of an AVX2 register.
	 */

	yi0 = _mm256_and_si256(ym38w,
		_mm256_setr_epi64x(
			dec64le(src),
			dec64le(src +  4) >> 6,
			dec64le(src +  9) >> 4,
			dec64le(src + 14) >> 2));
	yi1 = _mm256_and_si256(ym38w,
		_mm256_setr_epi64x(
			dec64le(src + 19),
			dec64le(src + 23) >> 6,
			dec64le(src + 28) >> 4,
			dec64le(src + 33) >> 2));

	/*
	 * Useful facts:
	 *
	 *   2^8 = 63 mod 193
	 *   2^16 = 109 mod 193
	 *   2^24 = 112 mod 193
	 *   2^32 = 108 mod 193
	 *   (x * 43465) >> 23 == floor(x / 193)  for all x in 0..61372.
	 *
	 * General processing is the following:
	 *
	 *  - Input is an integer z.
	 *  - We find a = z mod 193, which yields the next ciphertext value.
	 *  - z <- (z - a) / 193. This is an exact division, and the result
	 *    fits on 32 bits, so we can compute it by multiplying z
	 *    with 2425655105, which is the inverse of 193 modulo 2^32.
	 *    We can even do that by working with z mod 2^32.
	 */

	/*
	 * First element: split z into its five bytes, multiply them by
	 * the corresponding factor, and add.
	 */
	ya0 = _mm256_add_epi64(
		_mm256_and_si256(yi0, ym8w),
		_mm256_add_epi64(
			_mm256_add_epi64(
				_mm256_mullo_epi16(y63,
					_mm256_and_si256(_mm256_srli_epi64(
						yi0, 8), ym8w)),
				_mm256_mullo_epi16(y109,
					_mm256_and_si256(_mm256_srli_epi64(
						yi0, 16), ym8w))),
			_mm256_add_epi64(
				_mm256_mullo_epi16(y112,
					_mm256_and_si256(_mm256_srli_epi64(
						yi0, 24), ym8w)),
				_mm256_mullo_epi16(y108,
					_mm256_srli_epi64(yi0, 32)))));
	ya1 = _mm256_add_epi64(
		_mm256_and_si256(yi1, ym8w),
		_mm256_add_epi64(
			_mm256_add_epi64(
				_mm256_mullo_epi16(y63,
					_mm256_and_si256(_mm256_srli_epi64(
						yi1, 8), ym8w)),
				_mm256_mullo_epi16(y109,
					_mm256_and_si256(_mm256_srli_epi64(
						yi1, 16), ym8w))),
			_mm256_add_epi64(
				_mm256_mullo_epi16(y112,
					_mm256_and_si256(_mm256_srli_epi64(
						yi1, 24), ym8w)),
				_mm256_mullo_epi16(y108,
					_mm256_srli_epi64(yi1, 32)))));

	/*
	 * Continue computation of the first element.
	 *
	 *   input: x <= 79479
	 *   x <- (x mod 2^10) + 59*floor(x / 2^10)    (x <= 5565)
	 *   x <- x - 193*(floor((x * 43465) / 2^23)   (x <= 192)
	 */
	ya0 = _mm256_add_epi32(_mm256_and_si256(ya0, ym10),
		_mm256_mullo_epi16(_mm256_srli_epi32(ya0, 10), y59));
	ya1 = _mm256_add_epi32(_mm256_and_si256(ya1, ym10),
		_mm256_mullo_epi16(_mm256_srli_epi32(ya1, 10), y59));
	ya0 = _mm256_sub_epi32(ya0, _mm256_mullo_epi16(y193,
		_mm256_srli_epi32(_mm256_mulhi_epu16(ya0, y43465), 7)));
	ya1 = _mm256_sub_epi32(ya1, _mm256_mullo_epi16(y193,
		_mm256_srli_epi32(_mm256_mulhi_epu16(ya1, y43465), 7)));

	/*
	 * First elements produced, repack them into 32-bit elements in yd0:
	 *
	 *  yd0:  e0 e4 e1 e5 | e2 e6 e3 e7
	 */
	yd0 = _mm256_or_si256(ya0, _mm256_slli_epi64(ya1, 32));

	/*
	 * Subtract remainders and divide by 193. The division can be
	 * done with a multiplication by 2425655105 modulo 2^32. We also
	 * repack values in 32-bit elements (same order as yd0).
	 */
	yi0 = _mm256_blend_epi32(yi0, _mm256_slli_epi64(yi1, 32), 0xAA);
	yi0 = _mm256_mullo_epi32(_mm256_sub_epi32(yi0, yd0), y2425655105);

	/*
	 * Second iteration; elements go to yd1.
	 */
	ya0 = _mm256_add_epi32(
		_mm256_add_epi32(
			_mm256_and_si256(yi0, ym8),
			_mm256_mullo_epi16(y63,
				_mm256_and_si256(_mm256_srli_epi32(
					yi0, 8), ym8))),
		_mm256_add_epi32(
			_mm256_mullo_epi16(y109,
				_mm256_and_si256(_mm256_srli_epi32(
					yi0, 16), ym8)),
			_mm256_mullo_epi16(y112, _mm256_srli_epi32(yi0, 24))));
	yd1 = _mm256_sub_epi32(ya0, _mm256_mullo_epi16(y193,
		_mm256_srli_epi32(_mm256_mulhi_epu16(ya0, y43465), 7)));
	yi0 = _mm256_mullo_epi32(_mm256_sub_epi32(yi0, yd1), y2425655105);

	/*
	 * Third iteration; elements go to yd2. Input value is at most
	 * 7379470, which fits on 3 bytes.
	 */
	ya0 = _mm256_add_epi32(
		_mm256_add_epi32(
			_mm256_and_si256(yi0, ym8),
			_mm256_mullo_epi16(y63,
				_mm256_and_si256(_mm256_srli_epi32(
					yi0, 8), ym8))),
		_mm256_mullo_epi16(y109, _mm256_srli_epi32(yi0, 16)));
	yd2 = _mm256_sub_epi32(ya0, _mm256_mullo_epi16(y193,
		_mm256_srli_epi32(_mm256_mulhi_epu16(ya0, y43465), 7)));
	yi0 = _mm256_mullo_epi32(_mm256_sub_epi32(yi0, yd2), y2425655105);

	/*
	 * Fourth and fifth iterations: input value is at most
	 * floor((2^38-1) / 193^3) = 38235.
	 */
	yd4 = _mm256_srli_epi32(_mm256_mulhi_epu16(yi0, y43465), 7);
	yd3 = _mm256_sub_epi32(yi0, _mm256_mullo_epi16(y193, yd4));

	/*
	 * Result is correct if all elements of yd4 are lower than
	 * 193 at this point (by construction, the elements in yd0..yd3
	 * are all in the 0..192 range).
	 */
	yr = _mm256_cmpgt_epi32(y193, yd4);

	/*
	 * Repack all output elements into bytes.
	 *
	 * Low lanes:
	 *  yd0:   0  -  -  - 20  -  -  -  5  -  -  - 25  -  -  -
	 *  yd1:   1  -  -  - 21  -  -  -  6  -  -  - 26  -  -  -
	 *  yd2:   2  -  -  - 22  -  -  -  7  -  -  - 27  -  -  -
	 *  yd3:   3  -  -  - 23  -  -  -  8  -  -  - 28  -  -  -
	 *  yd4:   4  -  -  - 24  -  -  -  9  -  -  - 29  -  -  -
	 *
	 * High lanes:
	 *  yd0:  10  -  -  - 30  -  -  - 15  -  -  - 35  -  -  -
	 *  yd1:  11  -  -  - 31  -  -  - 16  -  -  - 36  -  -  -
	 *  yd2:  12  -  -  - 32  -  -  - 17  -  -  - 37  -  -  -
	 *  yd3:  13  -  -  - 33  -  -  - 18  -  -  - 38  -  -  -
	 *  yd4:  14  -  -  - 34  -  -  - 19  -  -  - 39  -  -  -
	 *
	 * We first merge yd0..yd3 together.
	 */
	yd0 = _mm256_or_si256(
		_mm256_or_si256(yd0, _mm256_slli_epi32(yd1, 8)),
		_mm256_or_si256(
			_mm256_slli_epi32(yd2, 16),
			_mm256_slli_epi32(yd3, 24)));

	/*
	 * Low lanes:
	 *  yd0:   0  1  2  3 20 21 22 23  5  6  7  8 25 26 27 28
	 *  yd4:   4  -  -  - 24  -  -  -  9  -  -  - 29  -  -  -
	 *
	 * High lanes:
	 *  yd0:  10 11 12 13 30 31 32 33 15 16 17 18 35 36 37 38
	 *  yd4:  14  -  -  - 34  -  -  - 19  -  -  - 39  -  -  -
	 */

	yd1 = _mm256_or_si256(
		_mm256_shuffle_epi8(yd0, _mm256_setr_epi8(
			0, 1, 2, 3, -1, 8, 9, 10,
			11, -1, -1, -1, -1, -1, -1, -1,
			9, 10, 11, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, 4, 5)),
		_mm256_shuffle_epi8(yd4, _mm256_setr_epi8(
			-1, -1, -1, -1, 0, -1, -1, -1,
			-1, 8, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, 8, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1)));
	yd2 = _mm256_or_si256(
		_mm256_shuffle_epi8(yd0, _mm256_setr_epi8(
			-1, -1, -1, -1, 4, 5, 6, 7,
			-1, 12, 13, 14, 15, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, 0, 1, 2, 3, -1, 8)),
		_mm256_shuffle_epi8(yd4, _mm256_setr_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1,
			4, -1, -1, -1, -1, 12, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, 0, -1)));
	yd3 = _mm256_or_si256(
		_mm256_shuffle_epi8(yd0, _mm256_setr_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1,
			6, 7, -1, 12, 13, 14, 15, -1,
			-1, -1, -1, -1, -1, -1, -1, -1)),
		_mm256_shuffle_epi8(yd4, _mm256_setr_epi8(
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, -1, -1, -1, -1, -1, -1,
			-1, -1, 4, -1, -1, -1, -1, 12,
			-1, -1, -1, -1, -1, -1, -1, -1)));

	/*
	 * Low lanes:
	 *  yd1:   0  1  2  3  4  5  6  7  8  9  -  -  -  -  -  -
	 *  yd2:   -  -  -  - 20 21 22 23 24 25 26 27 28 29  -  -
	 *  yd3:   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
	 *
	 * High lanes:
	 *  yd1:  16 17 18 19  -  -  -  -  -  -  -  -  -  - 30 31
	 *  yd2:   -  -  -  -  -  -  -  -  -  - 10 11 12 13 14 15
	 *  yd3:  32 33 34 35 36 37 38 39  -  -  -  -  -  -  -  -
	 */

	yd0 = _mm256_or_si256(yd1, _mm256_permute4x64_epi64(yd2, 0x4E));
	yd1 = _mm256_permute4x64_epi64(yd3, 0x4E);

	/*
	 * Output is now in yd0 (32 bytes) and yd3 (8 bytes). Values
	 * are still in the 0..192 range, which we must normalize to
	 * -96..+96.
	 */
	yd0 = _mm256_sub_epi8(yd0, y96n);
	yd1 = _mm256_sub_epi8(yd1, y96n);
	_mm256_storeu_si256((void *)(c +  0), yd0);
#if BAT_64
	enc64le(c + 32, _mm_cvtsi128_si64(_mm256_castsi256_si128(yd1)));
#else
	enc32le(c + 32, _mm_cvtsi128_si32(_mm256_castsi256_si128(yd1)));
	enc32le(c + 36, _mm_cvtsi128_si32(_mm256_castsi256_si128(
		_mm256_bsrli_epi128(yd1, 4))));
#endif

	return yr;
}

#endif

/* see inner.h */
size_t
bat_encode_ciphertext_769(void *out, size_t max_out_len,
	const int8_t *c, unsigned logn)
{
	size_t out_len, u, v, n;
	uint8_t *buf;
	uint32_t acc;
	int acc_len;

	n = (size_t)1 << logn;
	out_len = (((n / 5) * 38) + 8 * (n % 5) + 7) >> 3;
	if (out == NULL) {
		return out_len;
	}
	if (max_out_len < out_len) {
		return 0;
	}

	/*
	 * buf = output buffer
	 * v = current offset into output buffer
	 * acc = buffered bits
	 * acc_len = number of buffered bits (0..31)
	 * Note: unused bits in acc are zero.
	 */
	buf = out;
	v = 0;
	acc = 0;
	acc_len = 0;
	u = 0;

#if BAT_AVX2
	/*
	 * Encode sets of 8 groups of 5 values; each yields 38 bytes.
	 */
	for (; (u + 40) <= n; u += 40) {
		encode_ct_5x8(buf + v, c + u);
		v += 38;
	}
#endif

	/*
	 * Encode full groups of 5 values; each yields 38 bits.
	 */
	for (; (u + 5) <= n; u += 5) {
		uint32_t lo, hi;

		hi = encode_ct_5(&lo, c + u);
		acc |= (lo << acc_len);
		enc32le(buf + v, acc);
		v += 4;
		acc = (lo >> (31 - acc_len)) >> 1;
		acc |= hi << acc_len;
		acc_len += 6;
		if (acc_len >= 32) {
			enc32le(buf + v, acc);
			v += 4;
			acc_len -= 32;
			acc = hi >> (6 - acc_len);
		}
	}

	/*
	 * Encode remaining values, 8 bits per value.
	 */
	for (; u < n; u ++) {
		uint32_t w;

		w = (uint32_t)(c[u] + 96);
		acc |= w << acc_len;
		acc_len += 8;
		if (acc_len >= 32) {
			enc32le(buf + v, acc);
			v += 4;
			acc_len -= 32;
			acc = w >> (8 - acc_len);
		}
	}

	/*
	 * Flush buffered bits.
	 */
	while (acc_len > 0) {
		buf[v ++] = (uint8_t)acc;
		acc >>= 8;
		acc_len -= 8;
	}
	return v;
}

/* see inner.h */
size_t
bat_decode_ciphertext_769(int8_t *c, unsigned logn,
	const void *in, size_t max_in_len)
{
	const uint8_t *buf;
	uint32_t r, acc;
	size_t in_len, u, n;
	int acc_len;
#if BAT_AVX2
	__m256i yr;
#endif

	buf = in;
	n = (size_t)1 << logn;
	in_len = (((n / 5) * 38) + 8 * (n % 5) + 7) >> 3;
	if (max_in_len < in_len) {
		return 0;
	}
	u = 0;

#if BAT_AVX2
	yr = _mm256_set1_epi32(-1);
	for (; (u + 40) <= n; u += 40) {
		yr = _mm256_and_si256(yr, decode_ct_5x8(c + u, buf));
		buf += 38;
	}
	yr = _mm256_and_si256(yr, _mm256_bsrli_epi128(yr, 4));
	yr = _mm256_and_si256(yr, _mm256_bsrli_epi128(yr, 8));
	r = _mm_cvtsi128_si32(
		_mm_and_si128(
			_mm256_castsi256_si128(yr),
			_mm256_extracti128_si256(yr, 1))) & 1;
#else
	r = 1;
#endif

	/*
	 * Decode all blocks of 20 values (19 bytes exactly).
	 */
	for (; (u + 20) <= n; u += 20) {
		uint32_t w0, w1, w2, w3, w4;

		w0 = dec32le(buf);
		w1 = dec32le(buf + 4);
		w2 = dec32le(buf + 8);
		w3 = dec32le(buf + 12);
		w4 = dec24le(buf + 16);
		r &= decode_ct_5(c + u,
			w0, w1 & 0x3F);
		r &= decode_ct_5(c + u + 5,
			(w1 >> 6) | (w2 << 26), (w2 >> 6) & 0x3F);
		r &= decode_ct_5(c + u + 10,
			(w2 >> 12) | (w3 << 20), (w3 >> 12) & 0x3F);
		r &= decode_ct_5(c + u + 15,
			(w3 >> 18) | (w4 << 14), w4 >> 18);
		buf += 19;
	}

	/*
	 * Decode the full blocks of 5 values (38 bits). This requires
	 * some bit buffering management.
	 */
	acc = 0;
	acc_len = 0;
	for (; u + 5 <= n; u += 5) {
		uint32_t lo, hi;

		lo = acc;
		acc = dec32le(buf);
		buf += 4;
		lo |= (acc << acc_len);
		acc >>= acc_len;
		if (acc_len < 6) {
			acc |= (uint32_t)(*buf ++) << acc_len;
			acc_len += 8;
		}
		hi = acc & 0x3F;
		acc >>= 6;
		acc_len -= 6;
		r &= decode_ct_5(c + u, lo, hi);
	}

	/*
	 * Remaining values use 8 bits each.
	 */
	for (; u < n; u ++) {
		uint32_t w;

		if (acc_len < 8) {
			acc |= (uint32_t)(*buf ++) << acc_len;
		} else {
			acc_len -= 8;
		}
		w = acc & 0xFF;
		acc >>= 8;
		c[u] = (int)w - 96;
		r &= (uint32_t)(w - 193) >> 31;
	}

	/*
	 * There may be a few unused bits; if any is non-zero, then the
	 * input is invalid.
	 */
	r &= (uint32_t)(acc - 1) >> 31;

	return in_len & -(size_t)r;
}
