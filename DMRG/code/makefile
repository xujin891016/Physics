CC=g++
CFLAGS=-O2 -DNDEBUG 

EXE=test heisenberg

all: $(EXE)

test: ED.cpp RDM.cpp WR.cpp tensor.cpp test.cpp
	$(CC) -o test -I ./ ED.cpp RDM.cpp WR.cpp tensor.cpp test.cpp $(CFLAGS)

heisenberg: ED.cpp RDM.cpp WR.cpp tensor.cpp heisenberg.cpp
	$(CC) -o heisenberg -I ./ ED.cpp RDM.cpp WR.cpp tensor.cpp heisenberg.cpp $(CFLAGS)

clean:
	rm $(EXE) 2>&-
