CXX = g++

TARGET = MCContrails
SRCDIR = src
OBJDIR = obj

# Define compiler flags
# -Wall: Enable all standard warnings
# -Wextra: Enable extra warnings
# -g: Include debugging information
# -std=c++17: Use C++17 standard (adjust as needed, e.g., c++11, c++14, c++20)
CXXFLAGS = -Wall -Wextra -g -std=c++17 -I$(SRCDIR)
#CXXFLAGS = -Wall -Wextra -std=c++17 -I$(SRCDIR)

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

# Rule to clean up generated files
clean:
	@rm -rf $(OBJDIR) # Remove object files

# Phony targets: these are not actual files, but names for commands
.PHONY: all clean

