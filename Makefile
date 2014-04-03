default: binaries

Samtools-lib:
	cd samtools; make

bin:
	mkdir -p bin

Parallel:
	cd python/wiggletools; make

binaries: Parallel
	chmod 755 bin/*

test: tests

tests:
	cd test; python test.py

clean:
	cd samtools; make clean
	cd src; make clean
	rm bin/*
	rm lib/*
