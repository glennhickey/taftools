CWD:=$(shell pwd)
CXX ?= g++

CXXFLAGS := -O3 -Werror=return-type -std=c++14 -ggdb -g -MMD -MP $(CXXFLAGS)

TAF_DIR=$(CWD)/deps/taf

LIB_DEPS = $(TAF_DIR)/lib/stTaf.a $(TAF_DIR)/lib/sonLib.a
INC_DEPS = $(TAF_DIR)/inc/taf.h $(TAF_DIR)/inc/ond.h $(TAF_DIR)/inc/line_iterator.h
LIB_FLAGS = $(LIBS) $(LIB_DEPS)
INC_FLAGS = -I$(CWD) -I$(CWD)/deps/taf/inc -I$(CWD)/deps/taf/submodules/sonLib/lib

all: bin/taftools

bin/taftools: src/taftools.o src/uce.o $(LIB_DEPS)
	mkdir -p bin/
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o bin/taftools src/taftools.o src/uce.o $(LIB_FLAGS)

src/taftools.o: src/uce.hpp src/taftools.cpp $(INC_DEPS) 
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c src/taftools.cpp -o src/taftools.o $(INC_FLAGS)

src/uce.o: src/uce.hpp src/uce.cpp $(INC_DEPS)
	$(CXX) $(INCLUDE_FLAGS) $(CXXFLAGS) $(CPPFLAGS) -c src/uce.cpp -o src/uce.o $(INC_FLAGS)

clean:
	rm -rf bin/taftools src/*.o

static:
	CFLAGS="$${CFLAGS} -static" \
	CXXFLAGS="$${CXXFLAGS} -static" \
	${MAKE} all	
