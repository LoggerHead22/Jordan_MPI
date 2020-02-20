CC := mpicxx
CFLAGS := -std=c++14 -Wall -pg -O3 -ffast-math -march=native  -mfpmath=sse


OBJS := jordan.o

a.out: main.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lpthread
	
clean:
	rm -rf *.o *.out
	
%.o: %.cpp *.h
	$(CC) $(CFLAGS) -c $<
