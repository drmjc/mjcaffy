#!/bin/bash
#
# Determine the date that the CEL file was created, from the CEL file header
# eg "06/05/08 12:05:36"
#
# Mark Cowley, 2008-07-28
# 2009-11-12: modded to check for calvin/agcc files which I cannot yet read.
#

## V1.
# grep -m1 -a '^DatHeader' "$@" | egrep -o '[0-9]{2}/[0-9]{2}/[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}'

for arg in "$@"; do
	echo -n "$arg	"
	date=$(grep -m1 -a '^DatHeader' "$arg" | egrep -o '[0-9]{2}/[0-9]{2}/[0-9]{2}')
	if [ -z "$date" ]; then
		v5=$(grep -c -m1 -a 'affymetrix-calvin-intensity' "$arg")
		if [ $v5 -eq 1 ]; then
			echo "-	Can't get date/time from a Calvin CEL file"
		fi
	else
		echo $date
	fi
done

