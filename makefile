CC=g++
LD=g++

CFLAGS=-g -std=c++20 -O0
LDFLAGS=


ibm: main.o 
	$(LD) main.o $(LDFLAGS) -o ibm 

main.o: source/main.cpp 
	$(CC) source/main.cpp $(CFLAGS) -c

clean:
	rm *.o -f

