# Compiler-Einstellungen
WINDOWS_CC = x86_64-w64-mingw32-gcc
LINUX_CC = g++
CFLAGS = -Wall -Wextra -pedantic -std=c++11

# Verzeichnisse
LIB_DIR = lib_c
BIN_DIR = bin

# Dateien und Bibliothek
SOURCES = $(LIB_DIR)/bessel.cpp $(LIB_DIR)/resolve_q.cpp $(LIB_DIR)/resolves.cpp $(LIB_DIR)/resolve_srf_cs.cpp
SOURCES_ = bessel.cpp resolve_q.cpp resolves.cpp resolve_srf_cs.cpp
Windows_OBJECTS = $(SOURCES_:%.cpp=$(BIN_DIR)/Windows/%.o)
Linux_OBJECTS = $(SOURCES_:%.cpp=$(BIN_DIR)/Linux/%.o)
LIB_NAME = libcppcoil64

# Zielplattformen
TARGETS = Linux Windows

all: $(TARGETS)

Windows: $(BIN_DIR)/Windows/$(LIB_NAME).lib

Linux: $(BIN_DIR)/Linux/$(LIB_NAME).a

$(BIN_DIR)/Windows/$(LIB_NAME).lib: $(Windows_OBJECTS)
	ar -crs $@ $^

$(BIN_DIR)/Linux/$(LIB_NAME).a: $(Linux_OBJECTS)
	ar -crs $@ $^

$(BIN_DIR)/Windows/%.o: $(LIB_DIR)/%.cpp
	@mkdir -p $(BIN_DIR)/Windows
	$(WINDOWS_CC) $(CFLAGS) -c -o $@ $<

$(BIN_DIR)/Linux/%.o: $(LIB_DIR)/%.cpp
	@mkdir -p $(BIN_DIR)/Linux
	$(LINUX_CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf $(BIN_DIR)/*

.PHONY: all Windows Linux clean
