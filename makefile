# Makefile für das Kompilieren aller C++-Dateien im aktuellen Verzeichnis

# Compiler und Compiler-Flags
CXX := g++
CXXFLAGS := -Wall -Wextra -std=c++11

# Alle .cpp-Dateien im aktuellen Verzeichnis ermitteln
SRCS := $(wildcard *.cpp)

# Erzeugen der .o-Dateinamen aus den .cpp-Dateinamen
OBJS := $(patsubst %.cpp, %.o, $(SRCS))

# Name der ausführbaren Datei
TARGET := main

# Hauptziel (Standardziel)
all: $(TARGET)

# Regel zum Erzeugen des ausführbaren Programms
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Regel zum Kompilieren von .cpp zu .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Löscht alle generierten Dateien
clean:
	rm -f $(OBJS) $(TARGET)
