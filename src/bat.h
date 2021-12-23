#ifndef BAT_H__
#define BAT_H__

#include <stddef.h>
#include <stdint.h>

/*
 * For modulus qqq and degree nnn, the following types and functions are
 * defined:
 *
 *   bat_qqq_nnn_private_key
 *
 *      Private key structure; contains all private key elements, including
 *      a copy of the public key.
 *
 *   bat_qqq_nnn_public_key
 *
 *      Public key structure. Contains only the public key.
 *
 *   bat_qqq_nnn_ciphertext
 *
 *      Ciphertext structure. Contains the ciphertext polynomial and
 *      the FO tag.
 *
 *   int bat_qqq_nnn_keygen(
 *           bat_qqq_nnn_private_key *sk, void *tmp, size_t tmp_len);
 *
 *      Generate a new key pair. Returned value is 0 on success, a negative
 *      value on error. Buffer tmp (tmp_len bytes) should be large enough
 *      (see the BAT_qqq_nnn_TMP_KEYGEN macro).
 *
 *   void bat_qqq_nnn_get_public_key(
 *           bat_qqq_nnn_public_key *pk, const bat_qqq_nnn_private_key *sk);
 *
 *      Get a copy of the public key from the private key.
 *
 *   size_t bat_qqq_nnn_encode_private_key(
 *           void *out, size_t max_out_len,
 *           const bat_qqq_nnn_private_key *sk, int short_format);
 *
 *      Encode the private key into bytes. If short_format is zero, then
 *      the "long format" is used (encoding contains f, g, F, G, w, and the
 *      generation seed). If short_format is non-zero, then the "short
 *      format" is used (encoding contains only F and the seed). The short
 *      format is much smaller, but requires more CPU and temporary RAM
 *      when decoding.
 *
 *      If out is NULL, then max_out_len is ignored, and the function
 *      returns the size (in bytes) that the encoded private key would have.
 *      Otherwise, if the encoded private key would be longer than
 *      max_out_len, then the function returns 0 and encodes nothing.
 *      Otherwise, the encoded private key is written into out, and its
 *      size (in bytes) is returned.
 *
 *   size_t bat_qqq_nnn_decode_private_key(
 *           bat_qqq_nnn_private_key *sk,
 *           const void *in, size_t max_in_len,
 *           void *tmp, size_t tmp_len);
 *
 *      Decode the private key from bytes. If the incoming bytes are
 *      invalid, or relate to different set of parameters, or max_in_len
 *      is shorter than the private key size (i.e. it was truncated),
 *      then this function returns 0. Otherwise, it returns the actual
 *      size (in bytes) of the encoded private key (which is not greater
 *      than max_in_len, but may be lower than max_in_len).
 *
 *      If the encoded key uses the long format, then tmp and tmp_len
 *      are ignored. If the encoded key uses the short format, then
 *      tmp (of size tmp_len bytes) is used for temporary storage; in
 *      that case, if the buffer is too short, then the function fails
 *      and returns 0. The BAT_qqq_nnn_TMP_DECODE_PRIV macro evaluates to
 *      the required minimal size.
 *
 *   size_t bat_qqq_nnn_encode_public_key(
 *           void *out, size_t max_out_len,
 *           const bat_qqq_nnn_public_key *pk);
 *
 *      Encode the public key into bytes.
 *
 *      If out is NULL, then max_out_len is ignored, and the function
 *      returns the size (in bytes) that the encoded public key would have.
 *      Otherwise, if the encoded public key would be longer than
 *      max_out_len, then the function returns 0 and encodes nothing.
 *      Otherwise, the encoded public key is written into out, and its
 *      size (in bytes) is returned.
 *
 *   size_t bat_qqq_nnn_decode_public_key(
 *           bat_qqq_nnn_public_key *pk,
 *           const void *in, size_t max_in_len);
 *
 *      Decode the public key from bytes. If the incoming bytes are
 *      invalid, or relate to different set of parameters, or max_in_len
 *      is shorter than the public key size (i.e. it was truncated),
 *      then this function returns 0. Otherwise, it returns the actual
 *      size (in bytes) of the encoded public key (which is not greater
 *      than max_in_len, but may be lower than max_in_len).
 *
 *   size_t bat_qqq_nnn_encode_ciphertext(
 *           void *out, size_t max_out_len,
 *           const bat_qqq_nnn_ciphertext *ct);
 *
 *      Encode the ciphertext into bytes.
 *
 *      If out is NULL, then max_out_len is ignored, and the function
 *      returns the size (in bytes) that the encoded ciphertext would have.
 *      Otherwise, if the encoded ciphertext would be longer than
 *      max_out_len, then the function returns 0 and encodes nothing.
 *      Otherwise, the encoded ciphertext is written into out, and its
 *      size (in bytes) is returned.
 *
 *   size_t bat_qqq_nnn_decode_ciphertext(
 *           bat_qqq_nnn_ciphertext *ct,
 *           const void *in, size_t max_in_len);
 *
 *      Decode the ciphertext from bytes. If the incoming bytes are
 *      invalid, or relate to different set of parameters, or max_in_len
 *      is shorter than the ciphertext size (i.e. it was truncated),
 *      then this function returns 0. Otherwise, it returns the actual
 *      size (in bytes) of the encoded ciphertext (which is not greater
 *      than max_in_len, but may be lower than max_in_len).
 *
 *   int bat_qqq_nnn_encapsulate(
 *           void *secret, size_t secret_len,
 *           bat_qqq_nnn_ciphertext *ct,
 *           const bat_qqq_nnn_public_key *pk,
 *           void *tmp, size_t tmp_len);
 *
 *      Perform a key encpasulation with the provided public key. The
 *      resulting shared secret is written into secret[], while the
 *      ciphertext is written into *ct. The shared secret length is
 *      arbitrary (it internally comes from a BLAKE2-based KDF) but
 *      of course the sender and receiver should agree on the length to
 *      use, depending on what the secret is for.
 *
 *      On success, 0 is returned; on error, a negative error code is
 *      returned and the secret value is not produced. If provided
 *      temporary buffer (tmp, of size tmp_len bytes) is too small, then
 *      BAT_ERR_NOSPACE is returned (see BAT_qqq_nnn_TMP_ENCAPS).
 *
 *   int bat_qqq_nnn_decapsulate(
 *           void *secret, size_t secret_len,
 *           const bat_qqq_nnn_ciphertext *ct,
 *           const bat_qqq_nnn_private_key *sk,
 *           void *tmp, size_t tmp_len);
 *
 *      Perform a key decpasulation with the provided ciphertext and
 *      private key. The resulting shared secret is written into
 *      secret[]. The shared secret length is arbitrary (it internally
 *      comes from a BLAKE2-based KDF) but of course the sender and
 *      receiver should agree on the length to use, depending on what
 *      the secret is for.
 *
 *      On success, 0 is returned; on error, a negative error code is
 *      returned and the secret value is not produced. Such errors are
 *      reported only for local technical reasons unrelated to the
 *      received ciphertext; e.g. BAT_ERR_NOSPACE is returned if the
 *      tmp[] buffer (of size tmp_len bytes) is returned. By
 *      construction of the algorithm, invalid ciphertext values lead to
 *      a recovered shared secret which is deterministic from the
 *      ciphertext and private key, but unpredictable by third parties;
 *      in such cases, this function reports a success (0).
 */

#define BAT_MK(q, n, lvl_bytes, htype) \
typedef struct { \
	int8_t f[n]; \
	int8_t g[n]; \
	int8_t F[n]; \
	int8_t G[n]; \
	int32_t w[n]; \
	htype h[n]; \
	uint8_t rr[32]; \
	uint8_t seed[32]; \
} bat_ ## q ## _ ## n ## _private_key; \
typedef struct { \
	htype h[n]; \
} bat_ ## q ## _ ## n ## _public_key; \
typedef struct { \
	int8_t c[n]; \
	uint8_t c2[lvl_bytes]; \
} bat_ ## q ## _ ## n ## _ciphertext; \
int bat_ ## q ## _ ## n ## _keygen(bat_ ## q ## _ ## n ## _private_key *sk, \
	void *tmp, size_t tmp_len); \
void bat_ ## q ## _ ## n ## _get_public_key( \
	bat_ ## q ## _ ## n ## _public_key *pk, \
	const bat_ ## q ## _ ## n ## _private_key *sk); \
size_t bat_ ## q ## _ ## n ## _encode_private_key( \
	void *out, size_t max_out_len, \
	const bat_ ## q ## _ ## n ## _private_key *sk, int short_format); \
size_t bat_ ## q ## _ ## n ## _decode_private_key( \
	bat_ ## q ## _ ## n ## _private_key *sk, \
	const void *in, size_t max_in_len, \
	void *tmp, size_t tmp_len); \
size_t bat_ ## q ## _ ## n ## _encode_public_key( \
	void *out, size_t max_out_len, \
	const bat_ ## q ## _ ## n ## _public_key *pk); \
size_t bat_ ## q ## _ ## n ## _decode_public_key( \
	bat_ ## q ## _ ## n ## _public_key *pk, \
	const void *in, size_t max_in_len); \
size_t bat_ ## q ## _ ## n ## _encode_ciphertext( \
	void *out, size_t max_out_len, \
	const bat_ ## q ## _ ## n ## _ciphertext *ct); \
size_t bat_ ## q ## _ ## n ## _decode_ciphertext( \
	bat_ ## q ## _ ## n ## _ciphertext *ct, \
	const void *in, size_t max_in_len); \
int bat_ ## q ## _ ## n ## _encapsulate( \
	void *secret, size_t secret_len, \
	bat_ ## q ## _ ## n ## _ciphertext *ct, \
	const bat_ ## q ## _ ## n ## _public_key *pk, \
	void *tmp, size_t tmp_len); \
int bat_ ## q ## _ ## n ## _decapsulate( \
	void *secret, size_t secret_len, \
	const bat_ ## q ## _ ## n ## _ciphertext *ct, \
	const bat_ ## q ## _ ## n ## _private_key *sk, \
	void *tmp, size_t tmp_len);

BAT_MK(128, 256, 10, uint8_t)
BAT_MK(257, 512, 16, uint16_t)
BAT_MK(769, 1024, 32, uint16_t)

#undef BAT_MK

/*
 * Macros for temporary buffer sizes.
 *
 * Each length is in bytes and accounts for an extra 7 bytes for internal
 * alignment adjustment. If the provided buffera already ensures proper
 * generic alignment (e.g. it was obtained from malloc()), then the length
 * can be 7 bytes smaller.
 */
#define BAT_128_256_TMP_KEYGEN          6151
#define BAT_128_256_TMP_DECODE_PRIV     6151
#define BAT_128_256_TMP_ENCAPS           775
#define BAT_128_256_TMP_DECAPS          2055

#define BAT_257_512_TMP_KEYGEN         12295
#define BAT_257_512_TMP_DECODE_PRIV    12295
#define BAT_257_512_TMP_ENCAPS          2055
#define BAT_257_512_TMP_DECAPS          4103

#define BAT_769_1024_TMP_KEYGEN        24583
#define BAT_769_1024_TMP_DECODE_PRIV   24583
#define BAT_769_1024_TMP_ENCAPS         4103
#define BAT_769_1024_TMP_DECAPS         8199

/*
 * Error codes.
 */

/* Decapsulation failed. */
#define BAT_ERR_DECAPS_FAILED   -1

/* Provided object (key or ciphertext) uses a different set of parameters
   (modulus and/or degree) than expected by the called function. */
#define BAT_ERR_WRONG_PARAMS    -2

/* Provided object (key or ciphertext) is invalidly encoded. */
#define BAT_ERR_BAD_ENCODING    -3

/* Provided temporary space has insufficient length for the requested
   operation. */
#define BAT_ERR_NOSPACE         -4

/* Random seeding from operating system failed. */
#define BAT_ERR_RANDOM          -5

/*
 * Tag bytes. Each encoded public key, private key or ciphertext starts
 * with a tag byte that identifies the object type and parameters.
 * General format is (most-to-least significant order):
 *
 *    t t q q n n n n
 *
 * with:
 *
 *  - tt = 00 for a private key (long format), 01 for a private key (short
 *    format), 10 for a public key, 11 for a ciphertext.
 *  - qq = 00 for q = 128, 01 for q = 257, 10 for q = 769.
 *  - nnnn = log2(n) where n is the degree (power of 2, up to 1024).
 */
#define BAT_128_256_TAG_PRIVKEY_LONG     0x08
#define BAT_128_256_TAG_PRIVKEY_SHORT    0x48
#define BAT_128_256_TAG_PUBKEY           0x88
#define BAT_128_256_TAG_CIPHERTEXT       0xC8

#define BAT_257_512_TAG_PRIVKEY_LONG     0x19
#define BAT_257_512_TAG_PRIVKEY_SHORT    0x59
#define BAT_257_512_TAG_PUBKEY           0x99
#define BAT_257_512_TAG_CIPHERTEXT       0xD9

#define BAT_769_1024_TAG_PRIVKEY_LONG    0x2A
#define BAT_769_1024_TAG_PRIVKEY_SHORT   0x6A
#define BAT_769_1024_TAG_PUBKEY          0xAA
#define BAT_769_1024_TAG_CIPHERTEXT      0xEA

#endif
