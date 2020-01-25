GUROBIDIR = ../../../../gurobi900
INC = $(GUROBIDIR)/linux64/include/
LIB = $(GUROBIDIR)/linux64/lib/ 

rural: dsu.cpp rural.cpp
	g++ -c dsu.cpp
	g++ -m64 -g -o rural rural.cpp dsu.o -I $(INC) -L $(LIB) -lgurobi_c++ -lgurobi90 -lm
