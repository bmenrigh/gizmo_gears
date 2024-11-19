CC=gcc

CFLAGS=-std=gnu17 -Wall -Wextra -march=native -O2
#CFLAGS=-std=gnu17 -Wall -Wextra -march=native -O2 -g
#CFLAGS=-std=c99 -Wall -Wextra -march=native -O0 -frounding-math
LIBS=-lm -lquadmath -lpng

main: fast_gg_render

fast_gg_render: fast_gg_render.o
	$(CC) $(CFLAGS) -o fast_gg_render fast_gg_render.c $(LIBS)

fast_gg_render.o: fast_gg_render.c
	$(CC) $(CFLAGS) -c fast_gg_render.c

clean:
	rm -f fast_gg_render
	rm -f *.o
	rm -f *~
	rm -f \#*

