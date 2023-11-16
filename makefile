ALL = main.o UNEDF.o HFBTHO.o HFBTHO_solver.o
Target = main
CXX = g++
Flags = -llapack -llapacke -lblas -Wall

$(Target) : $(ALL)
	$(CXX) -o $(Target) $(ALL) $(Flags)


clean:
	rm *.o *.mod
