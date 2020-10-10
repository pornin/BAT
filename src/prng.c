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
#define SystemFunction036   NTAPI SystemFunction036
#include <NTSecAPI.h>
#undef SystemFunction036
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
	/*
	 * We could try to optimize this code with some caching of the
	 * file descriptor, but this raises extra difficulties (this is
	 * hard to make thread-safe without dabbling with a mutex). It
	 * is simpler to assume that any Unix-like platform for which it
	 * is worth optimizing performance will also have a recent
	 * enough OS to use getentropy() (possibly as a wrapper around
	 * getrandom()).
	 */
	{
		int f;

		f = open("/dev/urandom", O_RDONLY | O_CLOEXEC);
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
	/*
	 * Nominally, a "Win32" implementation should use CryptoAPI
	 * (CryptAcquireContext(), then CryptGenRandom()) but this is
	 * quite inefficient and error prone. Since Windows XP and
	 * Windows Server 2003, the RtlGenRandom() function (from
	 * advapi32.dll) offers a much direct road to the OS RNG.
	 */
	if (RtlGenRandom(seed, len)) {
		return 1;
	}
#endif
	return 0;
}
