/*
 * PRNG and interface to the system RNG.
 */

#include <assert.h>

#include "inner.h"

/*
 * Include relevant system header files. For Win32, this will also need
 * linking with advapi32.dll, which we trigger with an appropriate #pragma.
 */
#if BAT_RAND_GETENTROPY
#include <unistd.h>
#endif
#if BAT_RAND_URANDOM
#include <sys/types.h>
#if !BAT_RAND_GETENTROPY
#include <unistd.h>
#endif
#include <fcntl.h>
#include <errno.h>
#endif
#if BAT_RAND_WIN32
#include <windows.h>
#include <wincrypt.h>
#pragma comment(lib, "advapi32")
#endif

/* see inner.h */
int
bat_get_seed(void *seed, size_t len)
{
	(void)seed;
	if (len == 0) {
		return 1;
	}
#if BAT_RAND_GETENTROPY
	if (getentropy(seed, len) == 0) {
		return 1;
	}
#endif
#if BAT_RAND_URANDOM
	{
		int f;

		f = open("/dev/urandom", O_RDONLY);
		if (f >= 0) {
			while (len > 0) {
				ssize_t rlen;

				rlen = read(f, seed, len);
				if (rlen < 0) {
					if (errno == EINTR) {
						continue;
					}
					break;
				}
				seed = (uint8_t *)seed + rlen;
				len -= (size_t)rlen;
			}
			close(f);
			if (len == 0) {
				return 1;
			}
		}
	}
#endif
#if BAT_RAND_WIN32
	{
		HCRYPTPROV hp;

		if (CryptAcquireContext(&hp, 0, 0, PROV_RSA_FULL,
			CRYPT_VERIFYCONTEXT | CRYPT_SILENT))
		{
			BOOL r;

			r = CryptGenRandom(hp, (DWORD)len, seed);
			CryptReleaseContext(hp, 0);
			if (r) {
				return 1;
			}
		}
	}
#endif
	return 0;
}
