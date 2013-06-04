ALL: mrsfast

LDFLAGS=-s -static
LIBS=-L/net/gs/vol3/software/modules-sw/zlib/1.2.6/Linux/RHEL6/x86_64/lib/ -lm -lz
CFLAGS= -O2 -s 
DEBUG= -g

#mrsfast: baseFAST.o MrsFAST.o Common.o CommandLineParser.o RefGenome.o HashTable.o  Reads.o Output.o
#	gcc $^ -o $@ ${LIBS}
#	rm *.o

mrsfast: baseFAST.o MrsFAST.o Common.o CommandLineParser.o RefGenome.o HashTable.o  Reads.o Output.o
	gcc $^ -o $@ ${LDFLAGS} ${LIBS} ${CFLAGS}
	rm *.o
