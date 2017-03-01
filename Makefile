# Generic Makefile
# Author: Rixin Li

CC           = mpicc
CXX          = mpicxx

SRC_CXX      = $(wildcard *.cpp)
SOURCES      = $(SRC_CXX)
OBJECTS      = $(SRC_CXX:.cpp=.o)
INCL         = $(wildcard *.hpp)
EXEC         = plan

CXXFLAGS     = 
LDFLAGS      = 
LIBS         = -lmpi_cxx -lmpi
OPTIMIZE     = -O3 -m64 -Wall -std=c++11 -DMPI_ON -DNDEBUG

all: $(SOURCES) $(EXEC)

$(EXEC): $(OBJECTS)
	$(CXX) $(OPTIMIZE) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@
	rm $(OBJECTS)

%.o: %.cpp
	$(CXX) $(OPTIMIZE) $(CXXFLAGS) -c $< -o $@

.phony : clean

clean:
	rm -f $(OBJECTS) $(EXEC)
