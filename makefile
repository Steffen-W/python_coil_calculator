# Compiler-Einstellungen
CC = g++
CFLAGS = -Wall -Wextra -pedantic -std=c++11

# Verzeichnisse
LIB_DIR = lib_c
BIN_DIR = bin

# Dateien und Bibliothek
SOURCES = $(LIB_DIR)/bessel.cpp $(LIB_DIR)/resolve_q.cpp $(LIB_DIR)/resolves.cpp $(LIB_DIR)/resolve_srf_cs.cpp
OBJECTS = $(SOURCES:%.cpp=%.o)
LIB_NAME = libcppcoil64

all: $(BIN_DIR)/$(LIB_NAME).so

$(BIN_DIR)/$(LIB_NAME).so: $(OBJECTS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -shared -o $@ $^

%.o: %.cpp
	$(CC) $(CFLAGS) -fPIC -c -o $@ $<

clean:
	rm -rf $(BIN_DIR)/*.so $(LIB_DIR)/*.o

.PHONY: clean
