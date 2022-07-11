esercizio_2: esercizio_2.o random.o esercizio_2.hpp random.h
	g++ esercizio_2.o random.o -o esercizio_2
esercizio_2.o: esercizio_2.cpp random.cpp esercizio_2.hpp random.h
	g++ -c esercizio_2.cpp -o esercizio_2.o
random.o: random.cpp random.h
	g++ -c random.cpp -o random.o

clean:
	rm *.o

do:
	./esercizio_2


