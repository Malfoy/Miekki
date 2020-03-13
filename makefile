CC=g++
CFLAGS=-Wall -g -std=c++11 -pipe -funit-at-a-time -fopenmp -lz -Isparsepp -g
LDFLAGS=-lpthread -fopenmp -lz  -Isparsepp -g
OBJECTS=SIMDCompressionAndIntersection/bitpacking.o SIMDCompressionAndIntersection/integratedbitpacking.o
OBJECTS+=SIMDCompressionAndIntersection/simdbitpacking.o SIMDCompressionAndIntersection/usimdbitpacking.o
OBJECTS+=SIMDCompressionAndIntersection/simdintegratedbitpacking.o   SIMDCompressionAndIntersection/intersection.o
OBJECTS+=SIMDCompressionAndIntersection/varintdecode.o SIMDCompressionAndIntersection/streamvbyte.o
OBJECTS+=SIMDCompressionAndIntersection/simdpackedsearch.o SIMDCompressionAndIntersection/simdpackedselect.o
OBJECTS+=SIMDCompressionAndIntersection/frameofreference.o SIMDCompressionAndIntersection/for.o

EXEC=Miekki

DEBUG=0

ifeq ($(DEBUG), 1)
        CFLAGS+=-DDEBUG -Og
else
        CFLAGS+=-DNDEBUG -Ofast  -march=native -mtune=native
        LDFLAGS+=-flto
endif


all: $(EXEC)


Miekki: main.cpp miekki.o utils.o
	$(CC) -o $@ $^  $(LDFLAGS)  $(OBJECTS) -ISIMDCompressionAndIntersection/include

main.o: main.cpp
	$(CC) -o $@ -c $< $(CFLAGS) -ISIMDCompressionAndIntersection/include

miekki.o: Miekki.cpp utils.h
	$(CC) -o $@ -c $< $(CFLAGS) -ISIMDCompressionAndIntersection/include

utils.o: utils.cpp
	$(CC) -o $@ -c $< $(CFLAGS)



clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
