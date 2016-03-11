CXX	  = mpicxx

OPTIMIZE  = -O3 -m64 -Wall -g

SOURCES   = main.cpp global.cpp fileop.cpp tree.cpp
OBJECTS   = $(SOURCES:.cpp=.o)
INCL      = global.hpp fileop.hpp tree.hpp

CXXFLAGS  = -std=c++11 -I/opt/local/include/gcc5 -I/opt/local/include/libomp -I/opt/local/include/mpich-mp
LDFLAGS   = -L/opt/local/lib/gcc5 -L/opt/local/lib/libomp -L/opt/local/lib/mpich-mp
LIBS      = -lm -lmpicxx -lmpi -lpmpi

OPTIONS   = -DMPI_ON

EXEC      = plato

all: $(SOURCES) $(EXEC)

$(EXEC): $(OBJECTS) 
	$(CXX) $(OPTIMIZE) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@
	rm $(OBJECTS)
.cpp.o:
	$(CXX) $(OPTIMIZE) $(OPTIONS) $(CXXFLAGS) -c $< -o $@

.PHONY : clean

clean:
	 rm -f $(OBJECTS) $(EXEC)

#end
