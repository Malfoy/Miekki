CC=g++
CFLAGS=-Wall -g -std=c++11 -pipe -funit-at-a-time -fopenmp -lz -Isparsepp -g
LDFLAGS=-lpthread -fopenmp -lz  -Isparsepp -g

EXEC=Miekki

ASSERTS ?= $(DEBUG)
ifeq ($(DEBUG), 1)
        CFLAGS+=-DDEBUG -Og
else
        CFLAGS+=-DNDEBUG -Ofast -flto -march=native -mtune=native
        LDFLAGS+=-flto
endif


all: $(EXEC)



Miekki: main.o  miekki.o utils.o
	$(CC) -o $@ $^ $(LDFLAGS)

main.o: main.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

miekki.o: Miekki.cpp utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

utils.o: utils.cpp
	$(CC) -o $@ -c $< $(CFLAGS)



clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
