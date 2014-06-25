default: Parallel

Samtools-lib:
	cd samtools; make

Wiggletools: Samtools-lib
	cd src; make -e

Parallel: Wiggletools
	cd python/wiggletools; make

test: tests

tests:
	cd test; python test.py

clean:
	cd samtools; make clean
	cd src; make clean
	rm bin/*
	rm lib/*
