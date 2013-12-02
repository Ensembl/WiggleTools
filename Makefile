default: Samtools Wiggletools Parallel doc

Samtools:
	cd samtools; make

Wiggletools:
	cd src; make

Parallel:
	cd parallel; make

doc: Manual.pdf

Manual.pdf:
	cd doc; make

test: tests

tests:
	cd test; python test.py

clean:
	cd samtools; make clean
	cd src; make clean
	rm bin/*
	rm lib/*
	rm *.pdf
