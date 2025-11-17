CXX = g++
CXXFLAGS = -std=c++11 -O0 -g -fopenmp -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-maybe-uninitialized
INCLUDES = -I. -Ilibfileio -Ilibphasing
TARGET = ProgramPhasing

# Bibliothèques
LIBDIRS = libfileio libphasing
LIBS = libfileio/libfileio.a libphasing/libphasing.a

# Sources principales
MAIN_SOURCES = main.cpp phasing_program.cpp readinteger.cpp readreal.cpp readnegativereal.cpp
MAIN_OBJECTS = $(MAIN_SOURCES:.cpp=.o)

# Générateur de données
GENERATOR_TARGET = GenerateRandomData
GENERATOR_SOURCE = GenerateRandomData.cpp
GENERATOR_OBJECT = $(GENERATOR_SOURCE:.cpp=.o)

# Headers
HEADERS = types.h phasing_program.h libfileio/fileio.h libphasing/phasing.h

.PHONY: all clean libs generator

all: libs $(TARGET) $(GENERATOR_TARGET)

generator: $(GENERATOR_TARGET)

libs:
	@echo "Building libraries..."
	@for libdir in $(LIBDIRS); do \
		$(MAKE) -C $$libdir; \
	done

libfileio/libfileio.a:
	$(MAKE) -C libfileio

libphasing/libphasing.a:
	$(MAKE) -C libphasing

$(TARGET): $(MAIN_OBJECTS) $(LIBS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET) $(MAIN_OBJECTS) $(LIBS) -lm

$(GENERATOR_TARGET): $(GENERATOR_OBJECT)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(GENERATOR_TARGET) $(GENERATOR_OBJECT) -lm

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(MAIN_OBJECTS) $(TARGET) $(GENERATOR_OBJECT) $(GENERATOR_TARGET)
	@for libdir in $(LIBDIRS); do \
		$(MAKE) -C $$libdir clean; \
	done
