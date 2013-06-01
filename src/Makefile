SAMTOOLS=../../samtools/
CFLAGS=-g -Wall -O3
LIBDIR=../lib
BINDIR=../bin

default: clean lib binary
	
binary: lib wiggleTools.o
	mkdir -p ${BINDIR}
	${CC} ${CFLAGS} -L. -L${KENT_SRC}/lib/${MACHTYPE} -L${SAMTOOLS} -L${SAMTOOLS}/bcftools -L${LIBDIR} wiggleTools.o ${SAMTOOLS}/bam_plcmd.o ${SAMTOOLS}/sample.o ${SAMTOOLS}/bam2bcf.o ${SAMTOOLS}/errmod.o ${SAMTOOLS}/bam2bcf_indel.o -static -lwiggletools -ljkweb -lbam -lbcf -lm -lz -lpthread -lssl -lcrypto -ldl -o ${BINDIR}/wiggletools

lib: ${LIBDIR}/libwiggletools.a 

${LIBDIR}/libwiggletools.a: wiggleIterators.o wiggleReader.o bigWiggleReader.o wiggleMultiplexer.o wiggleReducers.o bedReader.o bigBedReader.o bamReader.o wiggleStreaming.o
	mkdir -p ${LIBDIR}
	ar rcs ${LIBDIR}/libwiggletools.a *.o

%.o: %.c; ${CC} ${CFLAGS} -I${KENT_SRC}/inc/ -I${SAMTOOLS} -D_PBGZF_USE -c $< -o $@

clean:
	rm -Rf *.o *.a wiggletools

samtools:
	cd ${SAMTOOLS}; make clean; make

