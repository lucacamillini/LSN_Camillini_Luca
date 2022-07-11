CC = mpicxx

main: main.o random.o ga.o main.hpp random.h ga.hpp
	$(CC) main.o random.o ga.o -o main
main.o: main.cpp random.cpp ga.cpp main.hpp random.h ga.hpp
	$(CC) -c main.cpp -o main.o

random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o
	
ga.o: ga.cpp ga.hpp random.h
	$(CC) -c ga.cpp -o ga.o

clean:
	rm *.o fitness.dat population.dat

do:
	rm *.dat
	mpiexec -np 4 ./main


