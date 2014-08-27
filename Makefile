CC=g++
CFLAGS=-lgsl -lgslcblas
DEPS = allele.h environment.h genotype.h locus.h modelFunctions.h population.h randomv.h stability.h

a:	runSimulation.cpp allele.cpp genotype.cpp locus.cpp environment.cpp modelFunctions.cpp population.cpp randomv.cpp stability.cpp
	$(CC) -O3 -std=c++0x -I eigen/ -o a.out runSimulation.cpp allele.cpp genotype.cpp locus.cpp environment.cpp modelFunctions.cpp population.cpp randomv.cpp stability.cpp $(CFLAGS)
debug:      runSimulation.cpp allele.cpp genotype.cpp locus.cpp environment.cpp modelFunctions.cpp population.cpp randomv.cpp stability.cpp
        $(CC) -g -O3 -std=c++0x -I eigen/ -o a.debug.out runSimulation.cpp allele.cpp genotype.cpp locus.cpp environment.cpp modelFunctions.cpp population.cpp randomv.cpp stability.cpp $(CFLAGS)
clean: 
	rm -f a.out
