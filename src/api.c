/*
 * This file is not meant to be compiled independently, but to be
 * included (with #include) by another C file. The caller shall
 * first define the Q, N and LOGN macros to relevant values (decimal
 * literal constants only).
 */

#if !defined Q || !defined N || !defined LOGN || !defined LVLBYTES
#error This module must not be compiled separately.
#endif

#include "bat.h"
#include "inner.h"

#define XCAT(x, y)    XCAT_(x, y)
#define XCAT_(x, y)   x ## y
#define XSTR(x)       XSTR_(x)
#define XSTR_(x)      #x

#define Zn(name)   XCAT(XCAT(XCAT(bat_, Q), XCAT(_, N)), XCAT(_, name))
#define ZN(name)   XCAT(XCAT(XCAT(BAT_, Q), XCAT(_, N)), XCAT(_, name))

/*
 * Ensure 8-byte alignment of the provided pointer. Returned value is the
 * aligned pointer. If, after alignment, the size is not at least equal
 * to min_tmp_len, then NULL is returned.
 */
static void *
tmp_align(void *tmp, size_t tmp_len, size_t min_tmp_len)
{
	unsigned off;

	if (tmp == NULL) {
		return NULL;
	}
	off = (8u - (unsigned)(uintptr_t)tmp) & 7u;
	if (tmp_len < off || (tmp_len - off) < min_tmp_len) {
		return NULL;
	}
	return (void *)((uintptr_t)tmp + off);
}

/*
 * Recompute the additional secret seed (rr) from the private key seed.
 */
static void
make_rr(Zn(private_key) *sk)
{
	shake_context sc;
	const char *label;

	label = "BAT-rr-" XSTR(Q) "-" XSTR(N) ":";
	shake_init(&sc, 256);
	shake_inject(&sc, label, strlen(label));
	shake_inject(&sc, sk->seed, sizeof sk->seed);
	shake_flip(&sc);
	shake_extract(&sc, sk->rr, sizeof sk->rr);
}

/*
 * Compute the hash function Hash_m(), used over the plaintext polynomial
 * 's' to generate the encryption seed. Output size matches the security
 * level.
 */
static void
hash_m(void *dst, const void *sbuf, size_t sbuf_len)
{
	shake_context sc;
	const char *label;

	label = "BAT-hash_m-" XSTR(Q) "-" XSTR(N) ":";
	shake_init(&sc, 256);
	shake_inject(&sc, label, strlen(label));
	shake_inject(&sc, sbuf, sbuf_len);
	shake_flip(&sc);
	shake_extract(&sc, dst, LVLBYTES);
}

/*
 * Compute the hash function Hash_s(), used over the message m to
 * generate the seed for encryption. Seed size is always exactly 32
 * bytes.
 */
static void
hash_s(void *seed, const void *m, size_t m_len)
{
	shake_context sc;
	const char *label;

	label = "BAT-hash_s-" XSTR(Q) "-" XSTR(N) ":";
	shake_init(&sc, 256);
	shake_inject(&sc, label, strlen(label));
	shake_inject(&sc, m, m_len);
	shake_flip(&sc);
	shake_extract(&sc, seed, 32);
}

/*
 * Sample the polynomial s using the provided seed (as part of
 * encryption).
 */
static void
sample_s(void *sbuf, size_t sbuf_len, const void *seed)
{
	shake_context sc;
	const char *label;

	label = "BAT-sample_s-" XSTR(Q) "-" XSTR(N) ":";
	shake_init(&sc, 256);
	shake_inject(&sc, label, strlen(label));
	shake_inject(&sc, seed, 32);
	shake_flip(&sc);
	shake_extract(&sc, sbuf, sbuf_len);
}

/*
 * Initialize a SHAKE256 context with the label for key derivation
 * and the plaintext s (in sbuf[]). Upon return, the context is in
 * production mode. The key derivation function is called Hash_k()
 * in the BAT specification.
 */
static void
init_kdf_good(shake_context *sc,
	const Zn(ciphertext) *ct, const void *m, size_t m_len)
{
	const char *label;

	label = "BAT-kdf-good-" XSTR(Q) "-" XSTR(N) ":";
	shake_init(sc, 256);
	shake_inject(sc, label, strlen(label));
	shake_inject(sc, m, m_len);
	shake_inject(sc, ct->c, sizeof ct->c);
	shake_inject(sc, ct->c2, sizeof ct->c2);
	shake_flip(sc);
}

/*
 * Initialize a SHAKE256 context with the label for alternate key
 * derivation, to be used on decapsulation failure. Upon return, the
 * context is in production mode. This function is called F() in the
 * BAT specification.
 */
static void
init_kdf_bad(shake_context *sc,
	const Zn(private_key) *sk, const Zn(ciphertext) *ct)
{
	const char *label;

	label = "BAT-kdf-bad-" XSTR(Q) "-" XSTR(N) ":";
	shake_init(sc, 256);
	shake_inject(sc, label, strlen(label));
	shake_inject(sc, sk->rr, sizeof sk->rr);
	shake_inject(sc, ct->c, sizeof ct->c);
	shake_inject(sc, ct->c2, sizeof ct->c2);
	shake_flip(sc);
}

/*
 * Make the secret value from the plaintext s.
 */
static void
make_secret(void *secret, size_t secret_len,
	const Zn(ciphertext) *ct, const void *m, size_t m_len)
{
	shake_context sc;

	init_kdf_good(&sc, ct, m, m_len);
	shake_extract(&sc, secret, secret_len);
}

/* see bat.h */
int
Zn(keygen)(Zn(private_key) *sk, void *tmp, size_t tmp_len)
{
	shake_context rng;
	uint8_t rng_seed[32];

	tmp = tmp_align(tmp, tmp_len, ZN(TMP_KEYGEN) - 7);
	if (tmp == NULL) {
		return BAT_ERR_NOSPACE;
	}
	if (!bat_get_seed(rng_seed, sizeof rng_seed)) {
		return BAT_ERR_RANDOM;
	}
	shake_init(&rng, 256);
	shake_inject(&rng, rng_seed, sizeof rng_seed);
	shake_flip(&rng);
	for (;;) {
		shake_extract(&rng, sk->seed, sizeof sk->seed);
		if (!bat_keygen_make_fg(sk->f, sk->g,
			(uint16_t *)sk->h, Q, LOGN,
			sk->seed, sizeof sk->seed, tmp))
		{
			continue;
		}
		if (!bat_keygen_solve_FG(sk->F, sk->G, sk->f, sk->g,
			Q, LOGN, tmp))
		{
			continue;
		}
		if (!bat_keygen_compute_w(sk->w,
			sk->f, sk->g, sk->F, sk->G, Q, LOGN, tmp))
		{
			continue;
		}
		make_rr(sk);
		return 0;
	}
}

/* see bat.h */
void
Zn(get_public_key)(Zn(public_key) *pk, const Zn(private_key) *sk)
{
	memmove(pk->h, sk->h, sizeof sk->h);
}

static size_t
get_privkey_length(const Zn(private_key) *sk, int short_format)
{
	if (short_format) {
		return 1 + sizeof(sk->seed) + bat_trim_i8_encode(
			NULL, 0, NULL, LOGN, bat_max_FG_bits[LOGN]);
	} else {
		return 1 + sizeof(sk->seed) + sizeof(sk->rr)
			+ bat_trim_i8_encode(NULL, 0,
				sk->f, LOGN, bat_max_fg_bits[LOGN])
			+ bat_trim_i8_encode(NULL, 0,
				sk->g, LOGN, bat_max_fg_bits[LOGN])
			+ bat_trim_i8_encode(NULL, 0,
				sk->F, LOGN, bat_max_FG_bits[LOGN])
			+ bat_trim_i8_encode(NULL, 0,
				sk->G, LOGN, bat_max_FG_bits[LOGN])
			+ bat_trim_i32_encode(NULL, 0,
				sk->w, LOGN, bat_max_w_bits[LOGN])
			+ XCAT(bat_encode_, Q)(NULL, 0, sk->h, LOGN);
	}
}

/* see bat.h */
size_t
Zn(encode_private_key)(void *out, size_t max_out_len,
	const Zn(private_key) *sk, int short_format)
{
	uint8_t *buf;
	size_t len, off, out_len;

	out_len = get_privkey_length(sk, short_format);
	if (out == NULL) {
		return out_len;
	}
	if (max_out_len < out_len) {
		return 0;
	}
	buf = out;
	if (short_format) {
		buf[0] = ZN(TAG_PRIVKEY_SHORT);
		memmove(&buf[1], sk->seed, sizeof sk->seed);
		off = 1 + sizeof sk->seed;
		len = bat_trim_i8_encode(buf + off, out_len - off,
			sk->F, LOGN, bat_max_FG_bits[LOGN]);
		if (len == 0) {
			/* This should never happen in practice. */
			return 0;
		}
		off += len;
		return off;
	} else {
		buf[0] = ZN(TAG_PRIVKEY_LONG);
		memmove(&buf[1], sk->seed, sizeof sk->seed);
		off = 1 + sizeof sk->seed;
		memmove(&buf[off], sk->rr, sizeof sk->rr);
		off += sizeof sk->rr;
		len = bat_trim_i8_encode(buf + off, out_len - off,
			sk->f, LOGN, bat_max_fg_bits[LOGN]);
		if (len == 0) {
			/* This should never happen in practice. */
			return 0;
		}
		off += len;
		len = bat_trim_i8_encode(buf + off, out_len - off,
			sk->g, LOGN, bat_max_fg_bits[LOGN]);
		if (len == 0) {
			/* This should never happen in practice. */
			return 0;
		}
		off += len;
		len = bat_trim_i8_encode(buf + off, out_len - off,
			sk->F, LOGN, bat_max_FG_bits[LOGN]);
		if (len == 0) {
			/* This should never happen in practice. */
			return 0;
		}
		off += len;
		len = bat_trim_i8_encode(buf + off, out_len - off,
			sk->G, LOGN, bat_max_FG_bits[LOGN]);
		if (len == 0) {
			/* This should never happen in practice. */
			return 0;
		}
		off += len;
		len = bat_trim_i32_encode(buf + off, out_len - off,
			sk->w, LOGN, bat_max_w_bits[LOGN]);
		if (len == 0) {
			/* This should never happen in practice. */
			return 0;
		}
		off += len;
		len = XCAT(bat_encode_, Q)(buf + off, out_len - off,
			sk->h, LOGN);
		if (len == 0) {
			/* This should never happen in practice. */
			return 0;
		}
		off += len;
		return off;
	}
}

/* see bat.h */
size_t
Zn(decode_private_key)(Zn(private_key) *sk, const void *in, size_t max_in_len,
	void *tmp, size_t tmp_len)
{
	const uint8_t *buf;
	size_t off, len;

	if (in == NULL || max_in_len == 0) {
		return 0;
	}
	buf = in;
	switch (buf[0]) {
	case ZN(TAG_PRIVKEY_SHORT):
		if (max_in_len < get_privkey_length(sk, 1)) {
			return 0;
		}
		memmove(sk->seed, buf + 1, sizeof sk->seed);
		off = 1 + sizeof sk->seed;
		len = bat_trim_i8_decode(sk->F, LOGN, bat_max_FG_bits[LOGN],
			buf + off, max_in_len - off);
		if (len == 0) {
			return 0;
		}
		off += len;
		tmp = tmp_align(tmp, tmp_len, ZN(TMP_DECODE_PRIV) - 7);
		if (tmp == NULL) {
			return 0;
		}
		if (!bat_keygen_make_fg(sk->f, sk->g,
			(uint16_t *)sk->h, Q, LOGN,
			sk->seed, sizeof sk->seed, tmp))
		{
			return 0;
		}
		if (!bat_keygen_rebuild_G(sk->G, sk->f, sk->g, sk->F,
			Q, LOGN, tmp))
		{
			return 0;
		}
		if (!bat_keygen_compute_w(sk->w,
			sk->f, sk->g, sk->F, sk->G, Q, LOGN, tmp))
		{
			return 0;
		}
		make_rr(sk);
		return off;
	case ZN(TAG_PRIVKEY_LONG):
		if (max_in_len < get_privkey_length(sk, 0)) {
			return 0;
		}
		memmove(sk->seed, buf + 1, sizeof sk->seed);
		off = 1 + sizeof sk->seed;
		memmove(sk->rr, buf + off, sizeof sk->rr);
		off += sizeof sk->rr;
		len = bat_trim_i8_decode(sk->f, LOGN, bat_max_fg_bits[LOGN],
			buf + off, max_in_len - off);
		if (len == 0) {
			return 0;
		}
		off += len;
		len = bat_trim_i8_decode(sk->g, LOGN, bat_max_fg_bits[LOGN],
			buf + off, max_in_len - off);
		if (len == 0) {
			return 0;
		}
		off += len;
		len = bat_trim_i8_decode(sk->F, LOGN, bat_max_FG_bits[LOGN],
			buf + off, max_in_len - off);
		if (len == 0) {
			return 0;
		}
		off += len;
		len = bat_trim_i8_decode(sk->G, LOGN, bat_max_FG_bits[LOGN],
			buf + off, max_in_len - off);
		if (len == 0) {
			return 0;
		}
		off += len;
		len = bat_trim_i32_decode(sk->w, LOGN, bat_max_w_bits[LOGN],
			buf + off, max_in_len - off);
		if (len == 0) {
			return 0;
		}
		off += len;
		len = XCAT(bat_decode_, Q)(sk->h, LOGN,
			buf + off, max_in_len - off);
		if (len == 0) {
			return 0;
		}
		off += len;
		return off;
	default:
		return 0;
	}
}

/* see bat.h */
size_t
Zn(encode_public_key)(void *out, size_t max_out_len, const Zn(public_key) *pk)
{
	uint8_t *buf;
	size_t out_len, len;

	out_len = 1 + XCAT(bat_encode_, Q)(NULL, 0, pk->h, LOGN);
	if (out == NULL) {
		return out_len;
	}
	if (max_out_len < out_len) {
		return 0;
	}
	buf = out;
	buf[0] = ZN(TAG_PUBKEY);
	len = XCAT(bat_encode_, Q)(buf + 1, max_out_len - 1, pk->h, LOGN);
	if (len == 0) {
		return 0;
	}
	return 1 + len;
}

/* see bat.h */
size_t
Zn(decode_public_key)(Zn(public_key) *pk, const void *in, size_t max_in_len)
{
	const uint8_t *buf;
	size_t len;

	if (max_in_len == 0) {
		return 0;
	}
	buf = in;
	if (buf[0] != ZN(TAG_PUBKEY)) {
		return 0;
	}
	len = XCAT(bat_decode_, Q)(pk->h, LOGN, buf + 1, max_in_len - 1);
	if (len == 0) {
		return 0;
	}
	return 1 + len;
}

/* see bat.h */
size_t
Zn(encode_ciphertext)(void *out, size_t max_out_len, const Zn(ciphertext) *ct)
{
	uint8_t *buf;
	size_t out_len, len, off;

	out_len = 1 + XCAT(bat_encode_ciphertext_, Q)(NULL, 0, ct->c, LOGN)
		+ sizeof ct->c2;
	if (out == NULL) {
		return out_len;
	}
	if (max_out_len < out_len) {
		return 0;
	}
	buf = out;
	buf[0] = ZN(TAG_CIPHERTEXT);
	off = 1;
	len = XCAT(bat_encode_ciphertext_, Q)(
		buf + off, max_out_len - off, ct->c, LOGN);
	if (len == 0) {
		return 0;
	}
	off += len;
	memcpy(buf + off, ct->c2, sizeof ct->c2);
	off += sizeof ct->c2;
	return off;
}

/* see bat.h */
size_t
Zn(decode_ciphertext)(Zn(ciphertext) *ct, const void *in, size_t max_in_len)
{
	const uint8_t *buf;
	size_t off, len;

	if (max_in_len < 1) {
		return 0;
	}
	buf = in;
	if (buf[0] != ZN(TAG_CIPHERTEXT)) {
		return 0;
	}
	off = 1;
	len = XCAT(bat_decode_ciphertext_, Q)(
		ct->c, LOGN, buf + off, max_in_len - off);
	if (len == 0) {
		return 0;
	}
	off += len;
	if (max_in_len - off < sizeof ct->c2) {
		return 0;
	}
	memcpy(ct->c2, buf + off, sizeof ct->c2);
	off += sizeof ct->c2;
	return off;
}

/* see bat.h */
int
Zn(encapsulate)(void *secret, size_t secret_len,
	Zn(ciphertext) *ct, const Zn(public_key) *pk,
	void *tmp, size_t tmp_len)
{
	tmp = tmp_align(tmp, tmp_len, ZN(TMP_ENCAPS) - 7);
	if (tmp == NULL) {
		return BAT_ERR_NOSPACE;
	}

	/*
	 * Encapsulation may theoretically fail if the resulting
	 * vector norm is higher than a specific bound. However, this
	 * is very rare (it cannot happen at all for q = 257). Thus,
	 * we expect not to have to loop. Correspondingly, it is more
	 * efficient to use the random seed from the OS directly.
	 */
	for (;;) {
		uint8_t m[LVLBYTES], seed[32], sbuf[SBUF_LEN(LOGN)];
		size_t u;

		/*
		 * Get a random message m from the OS.
		 */
		if (!bat_get_seed(m, sizeof m)) {
			return BAT_ERR_RANDOM;
		}

		/*
		 * Use Hash_s() to get the encryption seed.
		 */
		hash_s(seed, m, sizeof m);

		/*
		 * Derive the seed into the polynomial s.
		 */
		sample_s(sbuf, sizeof sbuf, seed);
#if N < 8
		/* For very reduced toy versions, we don't even have a
		   full byte, and we must clear the unused bits. */
		sbuf[0] &= (1u << N) - 1u;
#endif

		/*
		 * Compute c1. This may fail (rarely!) only for q = 769.
		 */
		if (!XCAT(bat_encrypt_, Q)(ct->c, sbuf, pk->h, LOGN, tmp)) {
			continue;
		}

		/*
		 * Make c2 = Hash_m(s) XOR m.
		 */
		hash_m(ct->c2, sbuf, sizeof sbuf);
		for (u = 0; u < sizeof m; u ++) {
			ct->c2[u] ^= m[u];
		}

		/*
		 * Produce the shared secret (output of a successful key
		 * exchange).
		 */
		make_secret(secret, secret_len, ct, m, sizeof m);

		return 0;
	}
}

/* see bat.h */
int
Zn(decapsulate)(void *secret, size_t secret_len,
	const Zn(ciphertext) *ct, const Zn(private_key) *sk,
	void *tmp, size_t tmp_len)
{
	uint8_t sbuf[SBUF_LEN(LOGN)], m[LVLBYTES];
	uint8_t sbuf_alt[SBUF_LEN(LOGN)], seed[32];
	int8_t *c_alt;
	size_t u;
	uint32_t d;
	uint8_t *secret_buf;
	shake_context sc;

	tmp = tmp_align(tmp, tmp_len, ZN(TMP_DECAPS) - 7);
	if (tmp == NULL) {
		return BAT_ERR_NOSPACE;
	}

	/*
	 * Inner decryption never fails (at least, it never reports
	 * a failure).
	 */
	XCAT(bat_decrypt_, Q)(sbuf, ct->c,
		sk->f, sk->g, sk->F, sk->G, sk->w, LOGN, tmp);

	/*
	 * From sbuf, we derive the mask that allows recovery of m
	 * out of the second ciphertext half (c2).
	 */
	hash_m(m, sbuf, sizeof sbuf);
	for (u = 0; u < sizeof m; u ++) {
		m[u] ^= ct->c2[u];
	}

	/*
	 * Decryption is valid if and only if we can reencrypt the
	 * obtained message m and get the exact same polynomial s
	 * and ciphertext c1.
	 */
	hash_s(seed, m, sizeof m);
	sample_s(sbuf_alt, sizeof sbuf_alt, seed);
#if N < 8
	sbuf_alt[0] &= (1u << N) - 1u;
#endif
	c_alt = tmp;
	tmp = tmp_align((void *)(c_alt + N), ZN(TMP_DECAPS) - N,
		ZN(TMP_ENCAPS) - 7);
	if (tmp == NULL) {
		/* This should never happen in practice. */
		return BAT_ERR_NOSPACE;
	}
	d = XCAT(bat_encrypt_, Q)(c_alt, sbuf_alt, sk->h, LOGN, tmp);
	d --;
	for (u = 0; u < sizeof sbuf; u ++) {
		d |= sbuf[u] ^ sbuf_alt[u];
	}
	for (u = 0; u < sizeof ct->c; u ++) {
		d |= (uint32_t)(ct->c[u] - c_alt[u]);
	}

	/*
	 * If encapsulation worked AND yielded the same ciphertext as
	 * received, then d == 0 at this point, and we want to produce
	 * the secret key as a hash of m. Otherwise, d != 0, and we
	 * must produce the secret as a hash of the received ciphertext
	 * and the secret value r (stored in sk->rr). We MUST NOT leak
	 * which was the case, and therefore we must always compute
	 * both hashes and perform constant-time conditional replacement.
	 */
	init_kdf_bad(&sc, sk, ct);
	shake_extract(&sc, secret, secret_len);
	init_kdf_good(&sc, ct, m, sizeof m);
	d = ((uint32_t)(d | -d) >> 31) - 1;
	secret_buf = secret;
	for (u = 0; u < secret_len; u ++) {
		uint8_t x;

		shake_extract(&sc, &x, 1);
		secret_buf[u] ^= d & (secret_buf[u] ^ x);
	}
	return 0;
}
