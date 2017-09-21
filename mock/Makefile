EXEC      = mock_lognormal
DIRS      = 

all: $(EXEC)

DIRS := # non standard library location e.g. /opt/local

# Add include and lib to DIRS
INCLDIRS  = $(foreach dir, $(DIRS), -I$(dir)/include)
LIBDIRS   = $(foreach dir, $(DIRS), -L$(dir)/lib) 
CXXFLAGS := $(INCLDIRS)

OBJS1 = mock_lognormal.o grid.o power.o cosmology.o power_spectrum.o
LIBS1 = $(LIBDIRS) -lboost_program_options -lfftw3 -lgsl -lgslcblas

mock_lognormal: $(OBJS1)
	$(CXX) $(OBJS1) $(LIBS1) -o $@

.PHONY: clean
clean :
	rm -f $(EXEC) $(OBJS1)
