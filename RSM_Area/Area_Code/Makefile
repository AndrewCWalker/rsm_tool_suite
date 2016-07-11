#CC = cc
CC = mpicc
CFLAGS = -g -pg -O2 -std=c99
LIBS = -lm -lgsl -lgslcblas

OBJS = area.o

area :  $(OBJS)
	$(CC) -o area $(OBJS) $(LIBS) $(CFLAGS)

area.o : area.c $(HEADS)
	$(CC) -c area.c $(CFLAGS)

clean :
	@rm -f *.o *~ area gmon.out
