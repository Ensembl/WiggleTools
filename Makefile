default: binaries

bin:
	mkdir -p bin

Wiggletools: bin
	cd src; make -e

Parallel: Wiggletools
	cd python/wiggletools; make

binaries: Parallel
	chmod 755 bin/*

test: tests

tests:
	cd test; python test.py

clean:
	cd src; make clean
	rm bin/*
	rm lib/*
