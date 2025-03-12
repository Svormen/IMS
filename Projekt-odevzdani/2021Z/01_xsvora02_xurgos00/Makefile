# IMS project
# Authors: Peter Urgos (xurgos00), Slavomir Svorada (xsvora02)
# 12/2021

CC=g++
CFLAGS=-g -std=c++17 -pedantic -Wall -Wextra -O2

MODULES_MAIN=main.o arburg.o

all: main

main: $(MODULES_MAIN)
	$(CC) $(CFLAGS) $^ -o main

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

run: main
	./main

.PHONY: clean
clean:
	rm -rf *.o main
