#!/bin/bash

cat ${*:3} | wigToBigWig - $1 $2
