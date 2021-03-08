#!/usr/bin/env bash
filename=$1
cat $filename | sed -E 's/(\-?\+?[0-9][0-9]) ([0-9][0-9]) ([0-9][0-9](\.[0-9][0-9]?)?)/\(\1, \2, \3\)/g' > out.txt

