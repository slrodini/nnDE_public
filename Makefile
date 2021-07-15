CC=gcc
CFLAGS=-Wall -O3

%.o: %.c
	$(CC) $(CFLAGS)  -c $^

main: llist.o main.o
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm *.o main
