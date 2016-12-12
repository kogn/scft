DEBUG = 0
MKL=1
DIM=2

INTEL_DIR= /usr/local/intel/
MATHEMATICA_DIR= /usr/local/Wolfram/Mathematica/10.3/
SRCDIR=./src/
OBJDIR=./obj/
CC = gcc 
CXX = g++
OPTS=-Ofast

LDFLAGS=
LIB= -L ../lib/
COMMON=  -I ../include/
CFLAGS=  -DNUM_THREADS=8 -DDIM=${DIM} -fopenmp


ifeq ($(DEBUG), 1)
	OPTS=-O0 -g
endif 

CFLAGS+=$(OPTS)

.PHONY: clean all

VPATH=./src/
EXEC = test_solver test_trans main test_config
all: obj data $(EXEC)

LDFLAGS += -llapacke -lsoft1

ifeq ($(MKL), 1)
	CFLAGS+= -DMKL
	COMMON += -I ${INTEL_DIR}/mkl/include/ -I ${INTEL_DIR}/mkl/include/fftw/
	LIB += -L ${INTEL_DIR}/interfaces/fftw3xc -L ${INTEL_DIR}/mkl/lib/intel64/
	LDFLAGS += -lfftw3xc_gnu -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
else
	LDFLAGS += -lfftw3_omp -lfftw3 
endif
COMMON += -I ${MATHEMATICA_DIR}/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/CompilerAdditions/
LIB += -L ${MATHEMATICA_DIR}/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64/CompilerAdditions/ \
	   -L${MATHEMATICA_DIR}/SystemFiles/Libraries/Linux-x86-64

LDFLAGS += -lopenblas  -lWSTP64i4 -lrt  -ldl -luuid  -lm

OBJS = $(patsubst %.c,%.o,$(addprefix $(OBJDIR), $(notdir $(wildcard $(SRCDIR)*.c))))
OBJS += $(patsubst %.cpp,%.o,$(addprefix $(OBJDIR), $(notdir $(wildcard $(SRCDIR)*.cpp))))
DEPS = $(wildcard $(SRCDIR)*.h $(SRCDIR)*.hpp) Makefile


test_config : test_config.cpp Config.cpp Config.h 
	$(CXX) $(filter-out %.h,$^) -o $@

main : $(filter-out $(addprefix $(OBJDIR),$(filter-out main.o,$(addsuffix .o,$(EXEC)))),$(OBJS))
	$(CXX) $(COMMON) $(CFLAGS) $^  $(LIB) $(LDFLAGS) -o $@

test_trans: $(filter-out $(addprefix $(OBJDIR),$(filter-out test_trans.o,$(addsuffix .o,$(EXEC)))),$(OBJS))
	$(CXX) $(COMMON) $(CFLAGS) $^ $(LIB) $(LDFLAGS) -o $@ 

test_solver: $(filter-out $(addprefix $(OBJDIR),$(filter-out test_solver.o,$(addsuffix .o,$(EXEC)))),$(OBJS))
	$(CXX) $(COMMON) $(CFLAGS) $^ $(LIB) $(LDFLAGS) -o $@ 

$(OBJDIR)%.o: %.c $(DEPS)
	$(CC) $(COMMON) $(CFLAGS) -c $< $(LIB)  -o $@ 

$(OBJDIR)%.o: %.cpp $(DEPS)
	$(CXX) $(COMMON) $(CFLAGS) -c $< $(LIB) -o $@

clean :
	rm $(OBJDIR)*.o
obj:
	mkdir -p obj
data:
	mkdir -p data 

