PROG=PC_ali

# Sources and objects
HEAD= allocate.h    	# header files
CSRC= main_PC_ali.c CV_statistics.c tree.c D_Cont.c clique_msa.c \
	nj.c align_mult.c print_pdb.c \
      Contact_divergence_aux.c PC_ali.c allocate.c align_ss.c \
      read_structures.c read_pdb_mammoth.c contact_list.c Print_pairwise.c \
      NeedlemanWunsch.c Profit_aux.c tm_score.c McLachlan_float.c
# seq_sim.c rank.c

COBJ= main_PC_ali.o CV_statistics.o tree.o D_Cont.o clique_msa.o \
      nj.o align_mult.o print_pdb.o \
      Contact_divergence_aux.o  PC_ali.o allocate.o align_ss.o \
      read_structures.o read_pdb_mammoth.o contact_list.o Print_pairwise.o \
      NeedlemanWunsch.o Profit_aux.o tm_score.o McLachlan_float.o
# seq_sim.o rank.o

####################
# Compiler flags
CC=gcc
FC=f77

FFLAG=-g -ffixed-line-length-132

#For debugging:
#CFLAGS= -g -std=c99 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
LDFLAGS=-lm
#CFLAGS= -g -std=c99 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -fsanitize=undefined
#LDFLAGS=-lm -fsanitize=undefined
#For production:
CFLAGS= -g -std=c99 -O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
#LDFLAGS=-lm 

#UNAME_S := $(shell uname -s)
#ARCH := $(shell arch)
# macos
#ifeq ($(UNAME_S), Darwin)
#      ifeq ($(ARCH), i386)
#            CFLAGS= -std=c99 -O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -march=core2
#      else
#            CFLAGS= -std=c99 -O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -mcpu=apple-m1
#      endif
#else
#      CFLAGS= -std=c99 -O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -march=core2
#endif
#For debugging & profilling
# CFLAGS= -Wall -std=c99 -pedantic -g -pg -fbounds-check -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
#LDFLAGS=-lm -L/usr/lib64/libg2c.so.0.0.0


all: $(FOBJ) $(COBJ)
	$(CC) $(FOBJ) $(COBJ) $(F77LIB) -o $(PROG) $(LDFLAGS)
	echo "Executable file generated:" $(PROG)

intel: CC=icc
intel: FC=ifort
intel: LDFLAGS=-lm
intel: CFLAGS=-Wall -static -O3 -xHOST -ipo -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -vec-report
intel: all

intel-static: CC=icc
intel-static: FC=ifort
intel-static: LDFLAGS=-lm -static
intel-static: CFLAGS=-Wall -static -O3 -xHOST -ipo -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -vec-report
intel-static: all

gnu: CC=g77
gnu: FC=f77
gnu: LDFLAGS=-lm -lg2c
gnu: CFLAGS=-Wall -O3 -march=nocona -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
gnu: FFLAGS= -O3 -ffixed-line-length-132
gnu: all

gnu-static: CC=g77
gnu-static: FC=f77
gnu-static: LDFLAGS=-lm -lg2c -static
gnu-static: CFLAGS=-Wall -O3 -march=nocona -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
gnu-static: FFLAGS= -O3 -ffixed-line-length-132
gnu-static: all

%.o: %.c
	$(CC) $(LDFLAGS) $(CFLAGS) -c $< -o $@

%.o: %.f
	$(FC)  $(FFLAG) -c $< -o $@

clean:
	rm -fr $(COBJ) $(FOBJ) $(PROG)
