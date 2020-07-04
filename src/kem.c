#include "inner.h"

/* see inner.h */
void
bat_seed_to_sbuf(uint8_t *sbuf, uint32_t q, unsigned logn,
	const void *seed, size_t seed_len)
{
	shake_context sc;
	uint8_t tt;
	unsigned nbits;
	const char *tag;

	shake_init(&sc, 256);
	switch (q) {
	case 128: tag = "BAT-enc-128:"; break;
	case 257: tag = "BAT-enc-257:"; break;
	default:  tag = "BAT-enc-769:"; break;
	}
	shake_inject(&sc, tag, 12);
	tt = logn;
	shake_inject(&sc, &tt, 1);
	shake_inject(&sc, seed, seed_len);
	shake_flip(&sc);
	nbits = 1u << logn;
	shake_extract(&sc, sbuf, (nbits + 7) >> 3);
	if (nbits < 8) {
		sbuf[0] &= (1 << nbits) - 1;
	}
}

/* see inner.h */
void
bat_sbuf_to_secret_and_mask(uint8_t *mask, void *secret, size_t secret_len,
	uint32_t q, unsigned logn, const uint8_t *sbuf)
{
	shake_context sc;
	uint8_t tt;
	const char *tag;

	shake_init(&sc, 256);
	switch (q) {
	case 128: tag = "BAT-kdf-128:"; break;
	case 257: tag = "BAT-kdf-257:"; break;
	default:  tag = "BAT-kdf-769:"; break;
	}
	shake_inject(&sc, tag, 12);
	tt = logn;
	shake_inject(&sc, &tt, 1);
	shake_inject(&sc, sbuf, SBUF_LEN(logn));
	shake_flip(&sc);
	shake_extract(&sc, mask, 32);
	shake_extract(&sc, secret, secret_len);
}
