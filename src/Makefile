#

CPC = g++
OPT = -std=c++11

CPP_FILES = Atom.cpp Options.cpp Parts.cpp Point.cpp main.cpp
BASENAME = $(basename $(CPP_FILES))
OBJS = $(addsuffix .o, $(BASENAME))

isparallel: $(OBJS)
	$(CPC) $(OPT) $(OBJS) -o $@

Point.o: Point.cpp Point.h
	$(CPC) $(OPT) -c $< -o $@

Atom.o: Atom.cpp Atom.h Point.h
	$(CPC) $(OPT) -c $< -o $@

Options.o: Options.cpp Options.h
	$(CPC) $(OPT) -c $< -o $@

Parts.o: Parts.cpp Parts.h Point.h Atom.h Options.h
	$(CPC) $(OPT) -c $< -o $@

main.o: main.cpp Point.h Atom.h Options.h Parts.h
	$(CPC) $(OPT) -c $< -o $@

.PHONY: clean
clean:
	rm $(OBJS)

