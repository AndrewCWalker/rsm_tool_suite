CC = cc
CFLAGS = -g -pg -O2 -std=c99
LIBS = -lm -lgsl -lgslcblas

OBJS = RSM.o

RSM :  $(OBJS)
	$(CC) -o RSM $(OBJS) $(LIBS) $(CFLAGS)

RSM.o: RSM.c $(HEADS)
	$(CC) -c RSM.c $(CFLAGS)