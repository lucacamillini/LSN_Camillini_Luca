esercizio1.1: esercizio1.1.o random.o esercizio1.1.h random.h
	g++ esercizio1.1.o random.o -o esercizio1.1
esercizio1.1.o: esercizio1.1.cpp random.cpp esercizio1.1.h random.h
	g++ -c esercizio1.1.cpp -o esercizio1.1.o

random.o: random.cpp random.h
	g++ -c random.cpp -o random.o

clean:
	rm *.o

do:
	./esercizio1.1


