#!/bin/bash

file=$1

sort -k1,1 -k2,2n -k3,3n -m ${file}x/* > $file
rm -Rf ${file}x
