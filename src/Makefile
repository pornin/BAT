.POSIX:

CC = clang
CFLAGS = -Wall -Wextra -Wshadow -Wundef -O2
LD = clang
LDFLAGS =
LIBS =

OBJ = api_128_256.o api_257_512.o api_769_1024.o codec.o fft.o fnr.o kem.o kem128.o kem257.o kem769.o keygen.o modqp.o prng.o sha3.o

all: test_bat speed_bat

clean:
	-rm -f $(OBJ) test_bat test_bat.o speed_bat speed_bat.o

test_bat: $(OBJ) test_bat.o
	$(LD) $(LDFLAGS) -o test_bat test_bat.o $(OBJ) $(LIBS)

speed_bat: speed_bat.o $(OBJ)
	$(LD) $(LDFLAGS) -o speed_bat speed_bat.o $(OBJ) $(LIBS)

api_128_256.o: api_128_256.c api.c bat.h inner.h sha3.h
	$(CC) $(CFLAGS) -c -o api_128_256.o api_128_256.c

api_257_512.o: api_257_512.c api.c bat.h inner.h sha3.h
	$(CC) $(CFLAGS) -c -o api_257_512.o api_257_512.c

api_769_1024.o: api_769_1024.c api.c bat.h inner.h sha3.h
	$(CC) $(CFLAGS) -c -o api_769_1024.o api_769_1024.c

codec.o: codec.c inner.h sha3.h
	$(CC) $(CFLAGS) -c -o codec.o codec.c

fft.o: fft.c inner.h sha3.h
	$(CC) $(CFLAGS) -c -o fft.o fft.c

fnr.o: fnr.c inner.h sha3.h
	$(CC) $(CFLAGS) -c -o fnr.o fnr.c

kem.o: kem.c inner.h sha3.h
	$(CC) $(CFLAGS) -c -o kem.o kem.c

kem128.o: kem128.c inner.h sha3.h
	$(CC) $(CFLAGS) -c -o kem128.o kem128.c

kem257.o: kem257.c modgen.c inner.h sha3.h
	$(CC) $(CFLAGS) -c -o kem257.o kem257.c

kem769.o: kem769.c modgen.c inner.h sha3.h
	$(CC) $(CFLAGS) -c -o kem769.o kem769.c

keygen.o: keygen.c inner.h sha3.h
	$(CC) $(CFLAGS) -c -o keygen.o keygen.c

modqp.o: modqp.c modgen.c inner.h sha3.h
	$(CC) $(CFLAGS) -c -o modqp.o modqp.c

prng.o: prng.c inner.h sha3.h
	$(CC) $(CFLAGS) -c -o prng.o prng.c

sha3.o: sha3.c inner.h sha3.h
	$(CC) $(CFLAGS) -c -o sha3.o sha3.c

speed_bat.o: speed_bat.c inner.h sha3.h
	$(CC) $(CFLAGS) -c -o speed_bat.o speed_bat.c

test_bat.o: test_bat.c inner.h sha3.h
	$(CC) $(CFLAGS) -c -o test_bat.o test_bat.c
