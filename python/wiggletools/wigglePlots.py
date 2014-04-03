import sys
import argparse
import subprocess
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import matplotlib.cm as cm

##############################################
## Multihistogram
##############################################

def make_histogram(infile, labels, out, format='pdf', normalised=True):
	counts = dict((X,[]) for X in labels)
	bin_positions = []

	file = open(infile)
	for line in file:
		items = line.strip().split('\t')
		bin_positions.append(float(items[0]))
		for i in range(len(labels)):
			counts[labels[i]].append(float(items[i+1]))
	file.close()

	assert len(labels) == len(counts.keys())

	values = numpy.array([bin_positions for X in labels]).T
	if normalised:
		weights = numpy.array([numpy.array(counts[X]) / sum(counts[X]) for X in labels]).T
	else:
		weights = numpy.array([counts[X] for X in labels]).T
	pyplot.hist(values, bins=len(bin_positions), weights=weights, label=labels)
	pyplot.legend()
	pyplot.savefig(out, format=format)

##############################################
## Profile curve
##############################################

def make_profile_curve(infile, out, format='pdf'):
	X = []
	Y = []
	file = open(infile)
	for line in file:
		items = line.strip().split('\t')
		X.append(float(items[0]))
		Y.append(float(items[1]))
	file.close()
	
	pyplot.plot(X, Y, '-')
	pyplot.savefig(out, format=format)

##############################################
## Profile matrix
##############################################

def make_profiles_matrix(infile, out, format='pdf'):
	M = []
	file = open(infile)
	for line in file:
		items = line.strip().split('\t')
		row = map(float, items[3:])
		M.append(row)
	file.close()

	pyplot.imshow(numpy.array(M),interpolation='nearest', cmap=cm.Greys_r)
	pyplot.savefig(out, format=format)

##############################################
## Overlap histogram
##############################################

def make_overlaps(infile, out, format='pdf'):
	values = dict()
	file = open(infile)
	for line in file:
		items = line.strip().split('\t')
		if len(items) != 4 and len(items) != 0:
			file.close()
			return
		if items[3] not in values:
			values[items[3]] = 0
		values[items[3]] += float(items[-1])
	file.close()

	keys = values.keys()
	labels = sorted(keys)
	ind = numpy.arange(len(keys))    # the x locations for the groups
	width = 0.35       # the width of the bars: can also be len(x) sequence
	pyplot.bar(ind, [values[X] for X in keys], width)
	pyplot.xticks(ind+width/2., keys)
	pyplot.savefig(out, format=format)

##############################################
## Convenience wrapper
##############################################

def make_plot(plot, infile, out, labels=None, format='pdf'):
	if plot == 'hist':
		make_histogram(infile, labels, out, format)
	elif plot == 'profile':
		make_profile_curve(infile, out, format)
	elif plot == 'profiles':
		make_profiles_matrix(infile, out, format)
	elif plot == 'overlaps':
		make_overlaps(infile, out, format)

##############################################
## Command line tool
##############################################

def get_options():
	parser = argparse.ArgumentParser(description='Wiggletools wrapper for plot generation.')
	parser.add_argument('--in', '-i', dest='infile', help='WiggleTools output',required=True)
	parser.add_argument('--labels',dest='labels',help='Labels to the different histograms', nargs='*')
	parser.add_argument('--out','-o', dest='out',help='Outfile', required=True)
	parser.add_argument('--plot','-p', dest='plot',help='Type of plot',choices=['hist','profile','profiles','overlaps'],required=True)
	options = parser.parse_args()
	return options

def main():
	options = get_options()
	make_plot(options.plot, options.infile, options.out, options.labels)

if __name__ == '__main__':
	main()
