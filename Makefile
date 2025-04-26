# Compiler and flags
CXX = c++
# CXXFLAGS = -O2 -std=c++17 -Wall -I./src/Healpix_3.83/include/healpix_cxx -I./src/cfitsio-4.6.2/include -L./src/Healpix_3.83/lib -L./src/cfitsio-4.6.2/lib -lhealpix_cxx -lsharp -lcfitsio

CXXFLAGS = -O2 -std=c++17 -Wall -I./src/Healpix_3.83/include/healpix_cxx -I./src/cfitsio-4.6.2/include
LDFLAGS = -L./src/Healpix_3.83/lib -L./src/cfitsio-4.6.2/lib -lhealpix_cxx -lsharp -lcfitsio -Wl,-rpath=./src/Healpix_3.83/lib -Wl,-rpath=./src/cfitsio-4.6.2/lib


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

# Compile rule for .cc -> .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJS) $(TARGET)
