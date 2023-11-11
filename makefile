ALL = main.o UNEDF.o HFBTHO.o HFBTHO_solver.o
Target = main
CXX = g++

$(Target) : $(ALL)
	$(CXX) -o $(Target) $(ALL)


clean:
	rm *.o *.mod
