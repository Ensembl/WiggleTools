CFLAGS=-g -Wall -O3 -std=gnu99
LIBDIR=../lib
BINDIR=../bin
TABIX_SRC=../../tabix
KENT_SRC=../../userApps/kent/src/
INC=-I${SAMTOOLS} -I${KENT_SRC}/inc -I${TABIX_SRC}
LIB_PATHS=-L${KENT_SRC}/lib/${MACHTYPE} -L${SAMTOOLS} -L${SAMTOOLS}/bcftools -L${LIBDIR} -L${TABIX_SRC}
LIBS= -lwiggletools ${KENT_SRC}/lib/local/jkweb.a -lbam -lbcf -ltabix -lz -lpthread -lssl -lcrypto -ldl -lgsl -lgslcblas -lm
OPTS=-D_PBGZF_USE
SAMTOOLS=../samtools

default: lib bin
	
bin: ${BINDIR}/wiggletools

${BINDIR}/wiggletools: ${LIBDIR}/libwiggletools.a wiggletools.o
	mkdir -p ${BINDIR}
	${CC} ${CFLAGS} ${LIB_PATHS} wiggletools.c ${LIBS} -o ${BINDIR}/wiggletools 

lib: ${LIBDIR}/libwiggletools.a 

${LIBDIR}/libwiggletools.a: wiggleIterator.o wigReader.o bigWiggleReader.o multiplexer.o reducers.o bedReader.o bigBedReader.o bamReader.o apply.o bigFileReader.o commandParser.o wigWriter.o statistics.o unaryOps.o multiSet.o setComparisons.o bufferedReader.o vcfReader.o bcfReader.o plots.o mWigWriter.o recycleBin.o fib.o
	mkdir -p ${LIBDIR}
	ar rcs ${LIBDIR}/libwiggletools.a *.o ${SAMTOOLS}/bam_plcmd.o ${SAMTOOLS}/sample.o ${SAMTOOLS}/bam2bcf.o ${SAMTOOLS}/errmod.o ${SAMTOOLS}/bam2bcf_indel.o

%.o: %.c; ${CC} ${CFLAGS} ${INC} ${OPTS} -c $< -o $@

clean:
	rm -Rf *.o *.a wiggletools
