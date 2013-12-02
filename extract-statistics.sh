#!/bin/bash


cat ./test5-0.vec | awk '{ if ($1 == "0") print $3"\t"$4 }' > all-collisions.tsv
cat ./test5-0.vec | awk '{ if ($1 == "1") print $3"\t"$4 }' > particle-collisions.tsv
cat ./test5-0.vec | awk '{ if ($1 == "2") print $3"\t"$4 }' > wall-collisions.tsv
cat ./test5-0.vec | awk '{ if ($1 == "3") print $3"\t"$4 }' > transfers.tsv
