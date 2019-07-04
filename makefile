CC=g++
CFLAGS= -Wall -Ofast -std=c++11  -flto -pipe -funit-at-a-time -fopenmp -lz -Isparsepp
LDFLAGS=-flto -lpthread -fopenmp -lz  -Isparsepp



EXEC=Miekki


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
