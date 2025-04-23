# Compiler and flags
CXX = c++
CXXFLAGS = -O2 -std=c++17 -Wall -Ihealpix_cxx -Icfitsio -Lhealpix_cxx -Lcfitsio -lhealpix_cxx -lcxxsupport  -lsharp -lcfitsio

# -I/Users/nathanielhayman/work/healpix_cuda/src/healpix_cxx/include -I/Users/nathanielhayman/work/healpix_cuda/src/cfitsio/include -L/Users/nathanielhayman/work/healpix_cuda/src/healpix_cxx/lib -L/Users/nathanielhayman/work/healpix_cuda/src/cfitsio/lib

# Source files
SRCS = register.cc

# Object files
OBJS = $(SRCS:.cc=.o)

# Output binary
TARGET = main

# Default target
all: $(TARGET)

# Link object files into executable
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET)

# Compile rule for .cc -> .o
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJS) $(TARGET)
