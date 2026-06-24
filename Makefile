# Makefile for C Statistics Library
# ==================================

# Compiler settings
CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O2
LDFLAGS = -lm

# Directories
SRCDIR = src
INCDIR = include
EXAMPLEDIR = examples
BUILDDIR = build

# Library sources and objects
LIB_SRCS = $(SRCDIR)/stats.c $(SRCDIR)/cstats.c
LIB_OBJS = $(BUILDDIR)/stats.o $(BUILDDIR)/cstats.o
LIB_NAME = libstats.a

# Example program
EXAMPLE_SRC = $(EXAMPLEDIR)/example_basic.c
EXAMPLE_BIN = example_basic

# Installation directories (can be overridden with PREFIX=/custom/path)
PREFIX ?= /usr/local
INSTALL_INCLUDE = $(PREFIX)/include
INSTALL_LIB = $(PREFIX)/lib

# Default target
.PHONY: all
all: $(BUILDDIR) $(LIB_NAME) $(EXAMPLE_BIN)

# Create build directory
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Compile stats.c
$(BUILDDIR)/stats.o: $(SRCDIR)/stats.c $(INCDIR)/stats.h
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $@

# Compile cstats.c
$(BUILDDIR)/cstats.o: $(SRCDIR)/cstats.c $(INCDIR)/cstats.h
	$(CC) $(CFLAGS) -I$(INCDIR) -c $< -o $@

# Create static library
$(LIB_NAME): $(LIB_OBJS)
	ar rcs $@ $^

# Build example program
$(EXAMPLE_BIN): $(EXAMPLE_SRC) $(LIB_NAME)
	$(CC) $(CFLAGS) -I$(INCDIR) $< -L. -lstats $(LDFLAGS) -o $@

# Run the example
.PHONY: example
example: $(EXAMPLE_BIN)
	./$(EXAMPLE_BIN)

# Install library and headers
.PHONY: install
install: $(LIB_NAME)
	install -d $(INSTALL_INCLUDE)
	install -m 644 $(INCDIR)/stats.h $(INSTALL_INCLUDE)/
	install -m 644 $(INCDIR)/cstats.h $(INSTALL_INCLUDE)/
	install -d $(INSTALL_LIB)
	install -m 644 $(LIB_NAME) $(INSTALL_LIB)/

# Clean up build artifacts
.PHONY: clean
clean:
	rm -rf $(BUILDDIR)
	rm -f $(LIB_NAME)
	rm -f $(EXAMPLE_BIN)
	rm -f *.o

# Help target
.PHONY: help
help:
	@echo "Available targets:"
	@echo "  all      - Build library and example (default)"
	@echo "  example  - Build and run the example program"
	@echo "  install  - Install headers and library to PREFIX (default: /usr/local)"
	@echo "  clean    - Remove all build artifacts"
	@echo "  help     - Show this help message"
	@echo ""
	@echo "Variables:"
	@echo "  PREFIX   - Installation prefix (default: /usr/local)"
	@echo "  CC       - C compiler (default: gcc)"
	@echo "  CFLAGS   - Compiler flags (default: -Wall -Wextra -std=c11 -O2)"
