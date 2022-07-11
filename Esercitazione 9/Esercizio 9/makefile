main: main.o random.o ga.o main.hpp random.h ga.hpp
	g++ main.o random.o ga.o -o main
main.o: main.cpp random.cpp ga.cpp main.hpp random.h ga.hpp
	g++ -c main.cpp -o main.o

random.o: random.cpp random.h
	g++ -c random.cpp -o random.o

ga.o: ga.cpp ga.hpp random.h
	g++ -c ga.cpp -o ga.o

clean:
	rm *.o fitness.dat population.dat

do:
	rm *.dat
	./main


