# Compiler and flags
CXX = c++
CXXFLAGS = -O2 -std=c++17 -Wall -Ihealpix_cxx -Icxx_helpers -Isharp -Icfitsio

# Source files
HEALPIX_SRCS = $(wildcard healpix_cxx/*.cc)
SRCS = register.cc $(HEALPIX_SRCS)

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
