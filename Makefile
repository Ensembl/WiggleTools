%.o: %.c; gcc -g -Wall -I${KENT_SRC}/inc/ -c $< -o $@

default: wiggleTools.o wiggleIterators.o wiggleReader.o bigWiggleReader.o wiggleMultiplexer.o wiggleReducers.o
	ar rcs libwiggletools.a *.o
	gcc -g -Wall -L. -L${KENT_SRC}/lib/x86_64/ -static -lwiggletools -ljkweb -lm -lz -lpthread -lssl -lcrypto -ldl -o wiggletools

clean:
	rm -Rf *.o *.a wiggletools

