#!/bin/bash
#
# Determine the chip type, from the CEL file header
# eg HuGene-1_0-st-v1
#
# Mark Cowley, 26/5/08
#
grep -m1 -a '^DatHeader' "$@" | grep -o '[a-zA-Z0-9_-]*\.1sq' | sed 's/\.1sq//'
