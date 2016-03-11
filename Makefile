CXX	  = mpicxx

OPTIMIZE  = -O3 -m64 -Wall -g

SOURCES   = readvtklis.cpp readvtk.cpp readlis.cpp fop.cpp global.cpp octree.cpp
OBJECTS   = $(SOURCES:.cpp=.o)
INCL      = readvtk.h readlis.h fop.h global.h octree.h

CXXFLAGS  = -std=c++11 -I/opt/local/include/gcc49 -I/opt/local/include/mpich-mp
LDFLAGS   = -L/opt/local/lib/gcc49 -L/opt/local/lib/mpich-mp
LIBS      = -lm

EXEC      = readvtklis

all: $(SOURCES) $(EXEC)

$(EXEC): $(OBJECTS) 
	$(CXX) $(OPTIMIZE) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@
	rm $(OBJECTS)
.cpp.o:
	$(CXX) $(OPTIMIZE) $(CXXFLAGS) -c $< -o $@

.PHONY : clean

clean:
	 rm -f $(OBJECTS) $(EXEC)

#end
