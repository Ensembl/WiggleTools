%.o: %.c; gcc -g -Wall -I${KENT_SRC}/inc/ -c $< -o $@

default: wiggleTools.o wiggleIterators.o wiggleReader.o bigWiggleReader.o
	gcc -g -Wall -L${KENT_SRC}/lib/x86_64/ *.o -static -ljkweb -lm -lz -lpthread -lssl -lcrypto -ldl -o wiggletools

clean:
	rm -Rf *.o wiggletools

