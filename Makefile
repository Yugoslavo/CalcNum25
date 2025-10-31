CC := gcc
# default: fast
COPT ?= -O3 -Wall
# enable with: make DEBUG=1
ifeq ($(DEBUG),1)
  COPT := -O0 -g -Wall -fno-omit-frame-pointer
endif

INCP := -I./primme/PRIMMESRC/COMMONSRC/
LIB  := -L./primme -lprimme -llapack -lblas -lm

default: main

main: main.o prob.o time.o interface_primme.o debug.o plotflux.o 
	$(CC) $(COPT) $^ -o $@ $(LIB)

main.o: main.c prob.h time.h interface_primme.h csr_io.h plotflux.h
	$(CC) $(COPT) -c $< -o $@ $(INCP)

%.o: %.c %.h
	$(CC) $(COPT) -c $< -o $@ $(INCP)

.PHONY: clean
clean:
	rm -f *.o main
