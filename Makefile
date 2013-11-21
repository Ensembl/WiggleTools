default: Samtools Wiggletools Bigwigs Parallel doc

Samtools:
	cd samtools; make

Wiggletools:
	cd src; make

Parallel:
	cd parallel; make

Bigwigs:
	cd bigwigs; make

doc: Manual.pdf

Manual.pdf:
	cd doc; make

clean:
	cd samtools; make clean
	cd bigwigs; make clean
	cd src; make clean
	rm bin/*
	rm lib/*
	rm *.pdf
