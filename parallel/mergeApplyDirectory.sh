#!/bin/bash

cat ${1}x/* | sort -k1,1 -k2,2n > $1
rm -Rf ${1}x
