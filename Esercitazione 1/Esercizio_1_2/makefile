esercizio1.2: esercizio1.2.o random.o esercizio1.2.hpp random.h
	g++ esercizio1.2.o random.o -o esercizio1.2
esercizio1.2.o: esercizio1.2.cpp random.cpp esercizio1.2.hpp random.h
	g++ -c esercizio1.2.cpp -o esercizio1.2.o

random.o: random.cpp random.h
	g++ -c random.cpp -o random.o

clean:
	rm *.o

do:
	./esercizio1.2


