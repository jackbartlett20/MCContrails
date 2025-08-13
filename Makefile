CXX = g++

TARGET = MCContrails
SRCDIR = src
OBJDIR = obj
INCLUDES = -I$(SRCDIR) -Iinclude

# Compiler flags explained
# -Wall: Enable all standard warnings
# -Wextra: Enable extra warnings
# -g: Include debugging information
# -std=c++17: Use C++17 standard (adjust as needed, e.g., c++11, c++14, c++20)
# -fopenmp: Enable OpenMP parallelisation
# -O0: No optimisation
# -O3: Aggresive optimisation

# Optimised option:
CXXFLAGS = -Wall -Wextra -std=c++20 -O3 $(INCLUDES) -fopenmp
# Debug option:
#CXXFLAGS = -Wall -Wextra -g -std=c++20 -O0 $(INCLUDES) -fopenmp

SRCS = $(wildcard $(SRCDIR)/*.cpp)

# Generate a list of object files from the source files
# Each .cpp file will correspond to a .o file in the object directory
OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))

# Default target: builds the executable
all: $(TARGET)

# Rule to link object files into the executable
$(TARGET): $(OBJS)
	@echo "Linking $(TARGET)..."
	$(CXX) $(OBJS) -o $@ $(CXXFLAGS)

# Generic rule to compile C++ source files into object files
# $@: the target of the rule (e.g., obj/main.o)
# $<: the first prerequisite (e.g., src/main.cpp)
# -c: compile only, do not link
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to clean up object files and executable
clean:
	@rm -rf $(OBJDIR) $(TARGET)

# Phony targets: these are not actual files, but names for commands
.PHONY: all clean

