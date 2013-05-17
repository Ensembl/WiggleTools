%.o: %.c; gcc -g -O3 -Wall -I${KENT_SRC}/inc/ -c $< -o $@

default: wiggleTools.o wiggleIterators.o wiggleReader.o bigWiggleReader.o
	gcc -g -O3 -Wall -L${KENT_SRC}/lib/x86_64/ *.o -static -ljkweb -lm -lz -lpthread -lssl -lcrypto -ldl -o wiggletools

clean:
	rm -Rf *.o wiggletools

