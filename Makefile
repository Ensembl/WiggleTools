default: Samtools Wiggletools Parallel

Samtools:
	cd samtools; make

Wiggletools:
	cd src; make

Parallel:
	cd parallel; make

test: tests

tests:
	cd test; python test.py

clean:
	cd samtools; make clean
	cd src; make clean
	rm bin/*
	rm lib/*
