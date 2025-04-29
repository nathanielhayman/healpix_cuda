# Compiler and flags
CXX = c++

CFITSIO_FOLDER = ./src/cfitsio
HEALPIX_FOLDER = ./src/Healpix_3.83

CXXFLAGS = -O2 -std=c++17 -Wall -I"$(HEALPIX_FOLDER)/include/healpix_cxx" -I"$(CFITSIO_FOLDER)/include"
LDFLAGS = -L"$(HEALPIX_FOLDER)/lib" -L"$(CFITSIO_FOLDER)/lib" -lhealpix_cxx -lsharp -lcfitsio

# Source files
SRCS = src/stb_impl.cpp src/register_naive.cpp 

# Object files
OBJDIR = build
OBJS = $(SRCS:src/%.cpp=$(OBJDIR)/%.o)

# Output binary
TARGET = main

# Default target
all: $(TARGET)

# Link object files into executable
$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(TARGET)

.PHONY: all clean

# Compile rule for .cpp -> .o
$(OBJDIR)/%.o: src/%.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -rf $(OBJDIR) $(TARGET)
