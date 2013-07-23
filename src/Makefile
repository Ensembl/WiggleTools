SAMTOOLS=../../samtools/
CFLAGS=-g -Wall -O3
LIBDIR=../lib
BINDIR=../bin
INC=-I${KENT_SRC}/inc/ -I${SAMTOOLS}
LIB_PATHS=-L${KENT_SRC}/lib/${MACHTYPE} -L${SAMTOOLS} -L${SAMTOOLS}/bcftools -L${LIBDIR}
LIBS=-lwiggletools -ljkweb -lbam -lbcf -lm -lz -lpthread -lssl -lcrypto -ldl
OPTS=-D_PBGZF_USE

default: lib bin
	
bin: ${BINDIR}/wiggletools

${BINDIR}/wiggletools: ${LIBDIR}/libwiggletools.a wiggleTools.o
	mkdir -p ${BINDIR}
	${CC} ${CFLAGS} ${LIB_PATHS} wiggleTools.c -static ${LIBS} -o ${BINDIR}/wiggletools 

lib: ${LIBDIR}/libwiggletools.a 

${LIBDIR}/libwiggletools.a: wiggleIterators.o wiggleReader.o bigWiggleReader.o wiggleMultiplexer.o wiggleReducers.o bedReader.o bigBedReader.o bamReader.o wiggleStreaming.o wiggleApply.o bigFileReader.o binaryFileReader.o commandParser.o
	mkdir -p ${LIBDIR}
	ar rcs ${LIBDIR}/libwiggletools.a *.o ${SAMTOOLS}/bam_plcmd.o ${SAMTOOLS}/sample.o ${SAMTOOLS}/bam2bcf.o ${SAMTOOLS}/errmod.o ${SAMTOOLS}/bam2bcf_indel.o

%.o: %.c; ${CC} ${CFLAGS} ${INC} ${OPTS} -c $< -o $@

clean:
	rm -Rf *.o *.a wiggletools
