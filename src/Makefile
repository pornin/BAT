.POSIX:

CC = clang

CFLAGS = -Wall -Wextra -Wshadow -Wundef -O3
# To allow clang to optimize for the local machine, add:
#    -march=native
# (this will activate automatic vectorization by clang, if the local machine
# can do it; but the binary may fail to run on older CPU)
#
# To use the AVX2-specific code, add:
#    -DBAT_AVX2
# (this will use all the code wil explicit AVX2 intrinsic functions; the
# binary won't run on a machine without AVX2)
#
# You can use -march=native and -DBAT_AVX2 at the same time, for the
# highest optimization:
# CFLAGS = -Wall -Wextra -Wshadow -Wundef -O3 -march=native -DBAT_AVX2

LD = clang
LDFLAGS =
LIBS =

OBJ = api_128_256.o api_257_512.o api_769_1024.o codec.o fft.o fnr.o kem128.o kem257.o kem769.o keygen.o modqp.o prng.o blake2b.o blake2s.o

all: test_bat speed_bat

clean:
	-rm -f $(OBJ) test_bat test_bat.o speed_bat speed_bat.o

test_bat: $(OBJ) test_bat.o
	$(LD) $(LDFLAGS) -o test_bat test_bat.o $(OBJ) $(LIBS)

speed_bat: speed_bat.o $(OBJ)
	$(LD) $(LDFLAGS) -o speed_bat speed_bat.o $(OBJ) $(LIBS)

api_128_256.o: api_128_256.c api.c bat.h inner.h blake2.h
	$(CC) $(CFLAGS) -c -o api_128_256.o api_128_256.c

api_257_512.o: api_257_512.c api.c bat.h inner.h blake2.h
	$(CC) $(CFLAGS) -c -o api_257_512.o api_257_512.c

api_769_1024.o: api_769_1024.c api.c bat.h inner.h blake2.h
	$(CC) $(CFLAGS) -c -o api_769_1024.o api_769_1024.c

codec.o: codec.c inner.h blake2.h
	$(CC) $(CFLAGS) -c -o codec.o codec.c

fft.o: fft.c inner.h blake2.h
	$(CC) $(CFLAGS) -c -o fft.o fft.c

fnr.o: fnr.c inner.h blake2.h
	$(CC) $(CFLAGS) -c -o fnr.o fnr.c

kem128.o: kem128.c inner.h blake2.h
	$(CC) $(CFLAGS) -c -o kem128.o kem128.c

kem257.o: kem257.c modgen.c inner.h blake2.h
	$(CC) $(CFLAGS) -c -o kem257.o kem257.c

kem769.o: kem769.c modgen.c inner.h blake2.h
	$(CC) $(CFLAGS) -c -o kem769.o kem769.c

keygen.o: keygen.c inner.h blake2.h
	$(CC) $(CFLAGS) -c -o keygen.o keygen.c

modqp.o: modqp.c modgen.c inner.h blake2.h
	$(CC) $(CFLAGS) -c -o modqp.o modqp.c

prng.o: prng.c inner.h blake2.h
	$(CC) $(CFLAGS) -c -o prng.o prng.c

blake2b.o: blake2b.c inner.h blake2.h
	$(CC) $(CFLAGS) -c -o blake2b.o blake2b.c

blake2s.o: blake2s.c inner.h blake2.h
	$(CC) $(CFLAGS) -c -o blake2s.o blake2s.c

speed_bat.o: speed_bat.c bat.h inner.h blake2.h
	$(CC) $(CFLAGS) -c -o speed_bat.o speed_bat.c

test_bat.o: test_bat.c bat.h inner.h blake2.h
	$(CC) $(CFLAGS) -c -o test_bat.o test_bat.c
