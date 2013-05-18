default: libwiggletools.a
	gcc -g -Wall -L. -L${KENT_SRC}/lib/x86_64/ -static -lwiggletools -ljkweb -lm -lz -lpthread -lssl -lcrypto -ldl -o wiggletools

libwiggletools.a: wiggleTools.o wiggleIterators.o wiggleReader.o bigWiggleReader.o wiggleMultiplexer.o wiggleReducers.o
	ar rcs libwiggletools.a *.o

%.o: %.c; gcc -g -Wall -I${KENT_SRC}/inc/ -c $< -o $@

clean:
	rm -Rf *.o *.a wiggletools

