# Compiler and flags
CXX = c++

CFITSIO_FOLDER = ./src/cfitsio
HEALPIX_FOLDER = ./src/Healpix_3.83

CXXFLAGS = -O2 -std=c++17 -Wall -I$(HEALPIX_FOLDER)/include/healpix_cxx -I$(CFITSIO_FOLDER)/include
LDFLAGS = -L$(HEALPIX_FOLDER)/lib -L$(CFITSIO_FOLDER)/lib -lhealpix_cxx -lsharp -lcfitsio -Wl,-rpath=$(HEALPIX_FOLDER)/lib -Wl,-rpath=$(CFITSIO_FOLDER)/lib


# -I/Users/nathanielhayman/work/healpix_cuda/src/healpix_cxx/include -I/Users/nathanielhayman/work/healpix_cuda/src/cfitsio/include -L/Users/nathanielhayman/work/healpix_cuda/src/healpix_cxx/lib -L/Users/nathanielhayman/work/healpix_cuda/src/cfitsio/lib

# Source files
SRCS = ./src/register_naive.cpp ./src/stb_impl.cpp

# Object files
OBJS = $(SRCS:.cc=.o)

# Output binary
TARGET = main

# Default target
all: $(TARGET)

# Link object files into executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $(TARGET)

# Compile rule for .cpp -> .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJS) $(TARGET)
