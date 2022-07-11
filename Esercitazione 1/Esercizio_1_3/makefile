CC = g++
CFLAGS = -g -Wall
TARGET = esercizio1.3

$(TARGET): $(TARGET).o random.o $(TARGET).hpp random.h
	g++ $(TARGET).o random.o -o $(TARGET)
$(TARGET).o: $(TARGET).cpp random.cpp $(TARGET).hpp random.h
	g++ -c $(TARGET).cpp -o $(TARGET).o

random.o: random.cpp random.h
	g++ -c random.cpp -o random.o

clean:
	rm *.o

do:
	./$(TARGET)


