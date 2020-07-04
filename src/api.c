/*
 * This file is not meant to be compiled independently, but to be
 * included (with #include) by another C file. The caller shall
 * first define the Q, N and LOGN macros to relevant values (decimal
 * literal constants only).
 */

#if !defined Q || !defined N || !defined LOGN
#error This module must not be compiled separately.
#endif

#include "bat.h"
#include "inner.h"

#define XCAT(x, y)    XCAT_(x, y)
#define XCAT_(x, y)   x ## y

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
		return 1 + sizeof(sk->seed)
			+ bat_trim_i8_encode(NULL, 0,
				sk->f, LOGN, bat_max_fg_bits[LOGN])
			+ bat_trim_i8_encode(NULL, 0,
				sk->g, LOGN, bat_max_fg_bits[LOGN])
			+ bat_trim_i8_encode(NULL, 0,
				sk->F, LOGN, bat_max_FG_bits[LOGN])
			+ bat_trim_i8_encode(NULL, 0,
				sk->G, LOGN, bat_max_FG_bits[LOGN])
			+ bat_trim_i16_encode(NULL, 0,
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
		len = bat_trim_i16_encode(buf + off, out_len - off,
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
		return off;
	case ZN(TAG_PRIVKEY_LONG):
		if (max_in_len < get_privkey_length(sk, 0)) {
			return 0;
		}
		memmove(sk->seed, buf + 1, sizeof sk->seed);
		off = 1 + sizeof sk->seed;
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
		len = bat_trim_i16_decode(sk->w, LOGN, bat_max_w_bits[LOGN],
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

	out_len = 1 + (sizeof ct->fo)
		+ XCAT(bat_encode_ciphertext_, Q)(NULL, 0, ct->c, LOGN);
	if (out == NULL) {
		return out_len;
	}
	if (max_out_len < out_len) {
		return 0;
	}
	buf = out;
	buf[0] = ZN(TAG_CIPHERTEXT);
	memmove(buf + 1, ct->fo, sizeof ct->fo);
	off = 1 + sizeof ct->fo;
	len = XCAT(bat_encode_ciphertext_, Q)(
		buf + off, max_out_len - off, ct->c, LOGN);
	if (len == 0) {
		return 0;
	}
	off += len;
	return off;
}

/* see bat.h */
size_t
Zn(decode_ciphertext)(Zn(ciphertext) *ct, const void *in, size_t max_in_len)
{
	const uint8_t *buf;
	size_t off, len;

	if (max_in_len < 1 + (sizeof ct->fo)) {
		return 0;
	}
	buf = in;
	if (buf[0] != ZN(TAG_CIPHERTEXT)) {
		return 0;
	}
	memmove(ct->fo, buf + 1, sizeof ct->fo);
	off = 1 + sizeof ct->fo;
	len = XCAT(bat_decode_ciphertext_, Q)(
		ct->c, LOGN, buf + off, max_in_len - off);
	if (len == 0) {
		return 0;
	}
	off += len;
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
		uint8_t seed[32], mask[32];
		uint8_t sbuf[SBUF_LEN(LOGN)];
		int i;

		if (!bat_get_seed(seed, sizeof seed)) {
			return BAT_ERR_RANDOM;
		}
		bat_seed_to_sbuf(sbuf, Q, LOGN, seed, sizeof seed);
		if (!XCAT(bat_encapsulate_, Q)(ct->c, sbuf, pk->h, LOGN, tmp)) {
			continue;
		}
		bat_sbuf_to_secret_and_mask(mask, secret, secret_len,
			Q, LOGN, sbuf);
		for (i = 0; i < 32; i ++) {
			ct->fo[i] = seed[i] ^ mask[i];
		}
		return 0;
	}
}

/* see bat.h */
int
Zn(decapsulate)(void *secret, size_t secret_len,
	const Zn(ciphertext) *ct, const Zn(private_key) *sk,
	void *tmp, size_t tmp_len)
{
	uint8_t sbuf[SBUF_LEN(LOGN)], sbuf2[SBUF_LEN(LOGN)], mask[32];
	size_t u;
	unsigned r;

	tmp = tmp_align(tmp, tmp_len, ZN(TMP_DECAPS) - 7);
	if (tmp == NULL) {
		return BAT_ERR_NOSPACE;
	}

	/*
	 * Decapsulation properly said never fails (at least, it never
	 * reports a failure).
	 */
	XCAT(bat_decapsulate_, Q)(sbuf, ct->c,
		sk->f, sk->g, sk->F, sk->G, sk->w, LOGN, tmp);

	/*
	 * Fujisaki-Okamoto transform: we use the mask derived from the
	 * obtained shared secret in order to decrypt the original seed,
	 * and recompute sbuf. If it does not yield the exact same sbuf
	 * as the one obtained from decapsulation, then this is an error.
	 * Comparison must be constant-time.
	 */
	bat_sbuf_to_secret_and_mask(mask, secret, secret_len, Q, LOGN, sbuf);
	for (u = 0; u < sizeof mask; u ++) {
		mask[u] ^= ct->fo[u];
	}
	bat_seed_to_sbuf(sbuf2, Q, LOGN, mask, sizeof mask);
	r = 0;
	for (u = 0; u < sizeof sbuf; u ++) {
		r |= sbuf[u] ^ sbuf2[u];
	}

	/*
	 * If there is an error, then we explicitly clear the derived key,
	 * to minimize risks that an application ignores the error and
	 * proceeds.
	 */
	if (r != 0) {
		memset(secret, 0, secret_len);
		return BAT_ERR_DECAPS_FAILED;
	} else {
		return 0;
	}
}
