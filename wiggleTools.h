#ifndef _WIGGLETOOLS_DEF_
#define _WIGGLETOOLS_DEF_

#include <fstream>

class WiggleIterator {
	protected:
		char * chrom;
		int start;
		int finish;
		double value;
		bool done;

	private:
		void print(FILE * );

	public:
		virtual void pop();
		bool current();
		WiggleIterator * operator + (WiggleIterator *);
		WiggleIterator * operator - (WiggleIterator *);
		WiggleIterator * operator * (WiggleIterator *);
		WiggleIterator * operator / (WiggleIterator *);
		WiggleIterator * operator ^ (double);
		WiggleIterator * x(double);
		char * getChrom();
		int getStart();
		int getFinish();
		double getValue();
		void setChrom(char *);
		void setStart(int);
		void setFinish(int);
		void setValue(double);
		void toFile(char *);
		void toStdout();
		double AUC();
};

class SumWiggleIterator : public WiggleIterator {
	WiggleIterator * iterA;
	WiggleIterator * iterB;

	public:
		void pop();
		SumWiggleIterator (WiggleIterator *, WiggleIterator*);
};

class ScaleWiggleIterator : public WiggleIterator {
	WiggleIterator * iter;
	float scalar;

	public:
		void pop();
		ScaleWiggleIterator (WiggleIterator *, double);
};

class LogWiggleIterator : public WiggleIterator {
	WiggleIterator * iter;
	float base;
	float baseLog;

	public:
		void pop();
		LogWiggleIterator (WiggleIterator *);
		LogWiggleIterator (WiggleIterator * , double);
};

class ExpWiggleIterator : public WiggleIterator {
	WiggleIterator * iter;
	float radix;
	float radixLog;

	public:
		void pop();
		ExpWiggleIterator (WiggleIterator *);
		ExpWiggleIterator (WiggleIterator *, double);
};

class PowerWiggleIterator: public WiggleIterator {
	WiggleIterator * iter;
	float exponant;

	public:
		void pop();
		PowerWiggleIterator (WiggleIterator *, double);
};

class ProductWiggleIterator : public WiggleIterator {
	WiggleIterator * iterA;
	WiggleIterator * iterB;

	public:
		void pop();
		ProductWiggleIterator (WiggleIterator * , WiggleIterator * );
};

class WiggleReader: public WiggleIterator {
	FILE * infile;
	size_t length;
	char * line;
	bool fixedStep;
	int step;
	int span;

	public:
		void pop();
		WiggleReader (char * );

	private:
		void readLine();
		void readHeader();
};
#endif
