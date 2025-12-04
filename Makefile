CC   := gcc
COPT ?= -O3 -Wall

# debugger
ifeq ($(DEBUG),1)
  COPT := -O0 -g -Wall -fno-omit-frame-pointer
endif

# --- PRIMME et ARPACK 

# Lieux ou se trouve Primme et arpack 
PRIMME_PATH := -I./primme/PRIMMESRC/COMMONSRC/
ARPACK_PATH := /usr/local

# include <arpack/arpack.h>
ARPACK_INC := -I$(ARPACK_PATH)/include

# --- JADAMILU et ILUPACK 

JADAMILU_DIR := ../JADAMILU
ARCHCOMP     := INT32GNU

# headers JADAMILU/ILUPACK
JADAMILU_INC := -I$(JADAMILU_DIR)/src/include

# répertoire des libs / objets JADAMILU
JADAMILU_LIBDIR := $(JADAMILU_DIR)/lib/$(ARCHCOMP)

# lib jadamilu + dépendances AMD/SuiteSparse
JADAMILU_LIBS := -L$(JADAMILU_LIBDIR) -ljadamilu -lamd -lsuitesparseconfig

# objets MC64 / MC21 compilés
JADAMILU_MC_OBJS := \
    $(JADAMILU_LIBDIR)/mc64s.o \
    $(JADAMILU_LIBDIR)/mc21s.o
#	

INCP := $(PRIMME_PATH) $(ARPACK_INC) $(JADAMILU_INC)

# Librairies
ARPACK_LIB := $(ARPACK_PATH)/lib/libarpack.a
BLASLAPACK := -llapack -lopenblas -lgfortran -lm
PRIMME_LIB := -L./primme -lprimme

# ordre important 
LIB := $(PRIMME_LIB) $(ARPACK_LIB) \
       $(JADAMILU_LIBS) $(JADAMILU_MC_OBJS) \
       $(BLASLAPACK) $(RPATH)

default: main

#Les objects
OBJS := main.o prob.o time.o interface_primme.o arpack_interface.o debug.o plotflux.o eulerprog.o 


main: $(OBJS)
	$(CC) $(COPT) $^ -o $@ $(LIB)

main.o: main.c prob.h time.h interface_primme.h arpack_interface.h csr_io.h plotflux.h eulerprog.h
	$(CC) $(COPT) -c $< -o $@ $(INCP)

%.o: %.c %.h
	$(CC) $(COPT) -c $< -o $@ $(INCP)

.PHONY: clean
clean:
	rm -f *.o main
