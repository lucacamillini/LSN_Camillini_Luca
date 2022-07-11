main: main.o random.o main.hpp random.h
	g++ main.o random.o -o main
main.o: main.cpp random.cpp main.hpp random.h
	g++ -c main.cpp -o main.o

random.o: random.cpp random.h
	g++ -c random.cpp -o random.o

clean:
	rm *.o

do:
	rm *.dat
	./main


