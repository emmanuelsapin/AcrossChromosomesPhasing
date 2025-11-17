CXX = g++
CXXFLAGS = -std=c++11 -O0 -g -fopenmp -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable
INCLUDES = -I..
TARGET = libfileio.a
SOURCES = fileio.cpp
OBJECTS = $(SOURCES:.cpp=.o)
HEADERS = fileio.h ../types.h ../readinteger.h ../readreal.h ../readnegativereal.h

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	ar rcs $(TARGET) $(OBJECTS)

# Explicit dependencies: recompile if source or headers change
fileio.o: fileio.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c fileio.cpp -o fileio.o

clean:
	rm -f $(OBJECTS) $(TARGET)

