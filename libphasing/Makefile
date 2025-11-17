CXX = g++
CXXFLAGS = -std=c++11 -O0 -g -fopenmp -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-maybe-uninitialized
INCLUDES = -I.. -I../libfileio
TARGET = libphasing.a
SOURCES = phasing.cpp
OBJECTS = $(SOURCES:.cpp=.o)
HEADERS = phasing.h ../types.h ../libfileio/fileio.h

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	ar rcs $(TARGET) $(OBJECTS)

# Explicit dependencies: recompile if source or headers change
# Force recompilation by making object depend on source and headers
phasing.o: phasing.cpp phasing.h ../types.h ../libfileio/fileio.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c phasing.cpp -o phasing.o

clean:
	rm -f $(OBJECTS) $(TARGET)

