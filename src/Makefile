CFLAGS=-g -Wall -O3 -std=gnu99
LIBDIR=../lib
BINDIR=../bin
INC=-I${HTSLIB_SRC}/htslib -I${KENT_SRC}/inc
LIB_PATHS=-L${KENT_SRC}/lib/${MACHTYPE} -L${LIBDIR} -L${HTSLIB_SRC}
LIBS= -lwiggletools ${KENT_SRC}/lib/local/jkweb.a ${HTSLIB_SRC}/libhts.a -lz -lpthread -lssl -lcrypto -ldl -lgsl -lgslcblas -lm
OPTS=-D_PBGZF_USE

default: lib bin
	
bin: ${BINDIR}/wiggletools

${BINDIR}/wiggletools: ${LIBDIR}/libwiggletools.a wiggletools.o
	mkdir -p ${BINDIR}
	${CC} ${CFLAGS} ${LIB_PATHS} wiggletools.c ${LIBS} -o ${BINDIR}/wiggletools 

lib: ${LIBDIR}/libwiggletools.a 

${LIBDIR}/libwiggletools.a: wiggleIterator.o wigReader.o bigWiggleReader.o multiplexer.o reducers.o bedReader.o bigBedReader.o bamReader.o apply.o bigFileReader.o commandParser.o wigWriter.o statistics.o unaryOps.o multiSet.o setComparisons.o bufferedReader.o vcfReader.o bcfReader.o plots.o mWigWriter.o recycleBin.o fib.o samReader.o
	mkdir -p ${LIBDIR}
	ar rcs ${LIBDIR}/libwiggletools.a *.o

%.o: %.c; ${CC} ${CFLAGS} ${INC} ${OPTS} -c $< -o $@

clean:
	rm -Rf *.o *.a wiggletools
