#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
using namespace std;

#include "wiggleTools.h"

static FILE * openOrFail(string filename, string description, string mode) {
	FILE * file;
	if (!(file = fopen(filename.c_str(), mode.c_str()))) {
		printf("Could not open %s %s, exiting...\n", (char *) description.c_str(), (char *) filename.c_str());
	}
	return file;
}

bool WiggleIterator::current() {
	return !done;
}

WiggleIterator * WiggleIterator::operator +(WiggleIterator * iterB) {
	return new SumWiggleIterator(this, iterB);
}

WiggleIterator * WiggleIterator::operator -(WiggleIterator * iterB) {
	return new SumWiggleIterator(this, iterB->x(-1));
}

WiggleIterator * WiggleIterator::operator *(WiggleIterator * iterB) {
	return new ProductWiggleIterator(this, iterB);
}

WiggleIterator * WiggleIterator::operator /(WiggleIterator * iterB) {
	return new ProductWiggleIterator(this, (*iterB) ^ -1.0);
}

WiggleIterator * WiggleIterator::operator ^(double scalar) {
	return new PowerWiggleIterator(this, scalar);
}

WiggleIterator * WiggleIterator::x(double scalar) {
	return new ScaleWiggleIterator(this, scalar);
}

WiggleIterator * sum(std::vector<WiggleIterator*> iters) {
	WiggleIterator * s = iters[0];
	for (int i = 1; i < iters.size(); i++)
		s = (*s) + iters[i];
	return s;
}

WiggleIterator * product(std::vector<WiggleIterator*> iters) {
	WiggleIterator * prod = iters[0];
	for (int i = 1; i < iters.size(); i++)
		prod = (*prod) * iters[i];
	return prod;
}

WiggleIterator * mean(std::vector<WiggleIterator*> iters) {
	return sum(iters)->x(1 / iters.size());
}

void WiggleIterator::print(FILE * out) {
	for ( ; !done; pop())
		fprintf(out, "%s\t%i\t%i\t%f\n", chrom, start, finish, value);
}

void WiggleIterator::toFile(char * filename) {
	FILE * file = openOrFail(filename, "output file", "w");
	print(file);
	fclose(file);
}

void WiggleIterator::toStdout() {
	print(stdout);
}


char * WiggleIterator::getChrom() {
	return chrom;
}

int WiggleIterator::getStart() {
	return start;
}

int WiggleIterator::getFinish() {
	return finish;
}

double WiggleIterator::getValue() {
	return value;
}	

void  WiggleIterator::setChrom(char * c) {
	chrom = c;
}
void  WiggleIterator::setStart(int i) {
	start = i;
}
 
void  WiggleIterator::setFinish(int i) {
	finish = i;
}

void  WiggleIterator::setValue(double d) {
	value = d;
}

void WiggleIterator::pop() {
	exit(1);
}

double WiggleIterator::AUC() {
	double total = 0;
	for(;!done; pop()) 
		total += (finish - start + 1) * value;
	return total;
}

//////////////////////////////////////////////////////
// Scaling operations
//////////////////////////////////////////////////////

ScaleWiggleIterator::ScaleWiggleIterator(WiggleIterator * i, double s) {
	iter = i;
	scalar = s;
	pop();
}

void ScaleWiggleIterator::pop() {
	if (iter->current()) {
		iter->pop();
		chrom = iter->getChrom();
		start = iter->getStart();
		finish = iter->getFinish();
		value = scalar * iter->getValue();
	} else {
		done = true;
	}
}


//////////////////////////////////////////////////////
// Log operations
//////////////////////////////////////////////////////

const double E = 2.71828128459045;

LogWiggleIterator::LogWiggleIterator(WiggleIterator * i) {
	iter = i;
	base = E;
	baseLog = 1;
	pop();
}

LogWiggleIterator::LogWiggleIterator(WiggleIterator * i, double s) {
	iter = i;
	base = s;
	baseLog = log(base);
	pop();
}

LogWiggleIterator log(WiggleIterator * iter) {
	return LogWiggleIterator(iter);
}

LogWiggleIterator log(WiggleIterator * iter, double base) {
	return LogWiggleIterator(iter, base);
}

void LogWiggleIterator::pop() {
	if (iter->current()) {
		iter->pop();
		chrom = iter->getChrom();
		start = iter->getStart();
		finish = iter->getFinish();
		value = log(iter->getValue()) / baseLog;
	} else {
		done = true;
	}
}

//////////////////////////////////////////////////////
// Exponentiation operations
//////////////////////////////////////////////////////

ExpWiggleIterator::ExpWiggleIterator(WiggleIterator * i, double s) {
	iter = i;
	radix = s;
	radixLog = log(radix);
	pop();
}

ExpWiggleIterator::ExpWiggleIterator(WiggleIterator * i) {
	iter = i;
	radix = E;
	radixLog = 1;
	pop();
}

ExpWiggleIterator exp(WiggleIterator * iter) {
	return ExpWiggleIterator(iter);
}

ExpWiggleIterator exp(WiggleIterator * iter, double radix) {
	return ExpWiggleIterator(iter, radix);
}

void ExpWiggleIterator::pop() {
	if (iter->current()) {
		iter->pop();
		chrom = iter->getChrom();
		start = iter->getStart();
		finish = iter->getFinish();
		value = exp(iter->getValue() * radixLog);
	} else {
		done = true;
	}
}

//////////////////////////////////////////////////////
// Power operations
//////////////////////////////////////////////////////

PowerWiggleIterator::PowerWiggleIterator(WiggleIterator * i, double s) {
	iter = i;
	exponant = s;
	pop();
}

void PowerWiggleIterator::pop() {
	if (iter->current()) {
		iter->pop();
		chrom = iter->getChrom();
		start = iter->getStart();
		finish = iter->getFinish();
		value = pow(iter->getValue(), exponant);
	} else {
		done = true;
	}
}

//////////////////////////////////////////////////////
// Sum operations
//////////////////////////////////////////////////////

SumWiggleIterator::SumWiggleIterator(WiggleIterator * a, WiggleIterator * b) {
	iterA = a;
	iterB = b;
	pop();
}

void SumWiggleIterator::pop() {
	if (!iterA->current() && !iterB->current()) {
		// All done
		done = true;
	} else if (!iterA->current()) {
		// A expired
		chrom = iterB->getChrom();
		start = iterB->getStart();
		finish = iterB->getFinish();
		value = iterB->getValue();
		iterB->pop();
	} else if (!iterB->current()) {
		// B expired
		chrom = iterA->getChrom();
		start = iterA->getStart();
		finish = iterA->getFinish();
		value = iterA->getValue();
		iterA->pop();
	} else {
		int chromDiff = strcmp(iterA->getChrom(), iterB->getChrom());

		if (chromDiff < 0) {
			chrom = iterA->getChrom();
			start = iterA->getStart();
			finish = iterA->getFinish();
			value = iterA->getValue();
			iterA->pop();
		} else if (chromDiff > 0) {
			chrom = iterB->getChrom();
			start = iterB->getStart();
			finish = iterB->getFinish();
			value = iterB->getValue();
			iterB->pop();
		} else {
			// Both iterators on the same chromosome:	
			chrom = iterA->getChrom();
			if (iterA->getStart() < iterB->getStart()) {
				start = iterA->getStart();	
				if (iterA->getFinish() < iterB->getStart()) {
					finish = iterA->getFinish();
					value = iterA->getValue();
					iterA->pop();
				} else {
					finish = iterB->getStart() - 1;
					value = iterA->getValue();
					iterA->setStart(iterB->getStart());
				}
			} else if (iterB->getStart() < iterA->getStart()) {
				start = iterB->getStart();	
				if (iterB->getFinish() < iterA->getStart()) {
					finish = iterB->getFinish();
					value = iterB->getValue();
					iterB->pop();
				} else {
					finish = iterA->getStart() - 1;
					value = iterB->getValue();
					iterB->setStart(iterA->getStart());
				}
			} else {
				start = iterA->getFinish();
				value = iterA->getValue() + iterB->getValue();
				if (iterA->getFinish() < iterB->getFinish()) {
					iterB->setStart(iterA->getFinish() + 1);
					iterA->pop();
				} else if (iterB->getFinish() < iterA->getFinish()) {
					iterA->setStart(iterB->getFinish() + 1);
					iterB->pop();
				} else {
					iterA->pop();
					iterB->pop();
				}
			}
		}
	}
}

//////////////////////////////////////////////////////
// Product operations
//////////////////////////////////////////////////////

ProductWiggleIterator::ProductWiggleIterator(WiggleIterator * a, WiggleIterator * b) {
	iterA = a;
	iterB = b;
	pop();
}

void ProductWiggleIterator::pop() {
	if (!iterA->current() || !iterB->current()) {
		done = true;
	} else {	
		int chromDiff = strcmp(iterA->getChrom(), iterB->getChrom());
		if (chromDiff < 0) {
			chrom = iterA->getChrom();
			start = iterA->getStart();
			finish = iterA->getFinish();
			value = 0;
			iterA->pop();
		} else if (chromDiff > 0) {
			chrom = iterB->getChrom();
			start = iterB->getStart();
			finish = iterB->getFinish();
			value = 0;
			iterB->pop();
		} else {
			// Both iterators on the same chromosome:	
			chrom = iterA->getChrom();
			if (iterA->getStart() < iterB->getStart()) {
				start = iterA->getStart();	
				if (iterA->getFinish() < iterB->getStart()) {
					finish = iterA->getFinish();
					value = 0;
					iterA->pop();
				} else {
					finish = iterB->getStart() - 1;
					value = 0;
					iterA->setStart(iterB->getStart());
				}
			} else if (iterB->getStart() < iterA->getStart()) {
				start = iterB->getStart();	
				if (iterB->getFinish() < iterA->getStart()) {
					finish = iterB->getFinish();
					value = 0;
					iterB->pop();
				} else {
					finish = iterA->getStart() - 1;
					value = 0;
					iterB->setStart(iterA->getStart());
				}
			} else {
				start = iterA->getFinish();
				value = iterA->getValue() * iterB->getValue();
				if (iterA->getFinish() < iterB->getFinish()) {
					iterB->setStart(iterA->getFinish() + 1);
					iterA->pop();
				} else if (iterB->getFinish() < iterA->getFinish()) {
					iterA->setStart(iterB->getFinish() + 1);
					iterB->pop();
				} else {
					iterA->pop();
					iterB->pop();
				}
			}
		}
	}
}

//////////////////////////////////////////////////////
// File Reader
//////////////////////////////////////////////////////

WiggleReader::WiggleReader(char * f) {
	if (strcmp(f, "-")) {
		infile = openOrFail(f, "input file", "r");
	} else {
		infile = stdin;
	}
	size_t length = 5000;
	char * line = (char *) calloc(length, 1);
	pop();
}	

void WiggleReader::pop() {
	if (getline(&line, &length, infile) != EOF) {
		sscanf(line, "%s\t%i\t%i\t%lf", chrom, &start, &finish, &value);
	} else {
		done = true;
		fclose(infile);
	}
}

//////////////////////////////////////////////////////
// Command line executable
//////////////////////////////////////////////////////

static void printHelp() {
	puts("Help!");
}

int main(int argc, char ** argv) {
	if (argc < 2 || strcmp(argv[1], "help") == 0) {
		printHelp();
		return 0;
	} else if (strcmp(argv[1], "add") == 0) {
		SumWiggleIterator(new WiggleReader(argv[2]), new WiggleReader(argv[3])).toStdout();
	} else if (strcmp(argv[1], "scale") == 0) {
		ScaleWiggleIterator(new WiggleReader(argv[2]), atoi(argv[3])).toStdout();
	}

	return 1;
}
