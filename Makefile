default: libwiggletools.a wiggleTools.c
	gcc -g -O3 -Wall -L. -L${KENT_SRC}/lib/x86_64/ wiggleTools.c -static -lwiggletools -ljkweb -lm -lz -lpthread -lssl -lcrypto -ldl -o wiggletools

libwiggletools.a: wiggleIterators.o wiggleReader.o bigWiggleReader.o wiggleMultiplexer.o wiggleReducers.o
	ar rcs libwiggletools.a *.o

%.o: %.c; gcc -g -O3 -Wall -I${KENT_SRC}/inc/ -c $< -o $@

clean:
	rm -Rf *.o *.a wiggletools

