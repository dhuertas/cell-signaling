#!/bin/bash

c=1

for i in {5..355..7}
do

	> diffusion-results/collisions-$c.tsv

	> diffusion-results/collision-position-$c-x.tsv
	> diffusion-results/collision-position-$c-y.tsv
	> diffusion-results/collision-position-$c-z.tsv

	> diffusion-results/position-$c-x.tsv
	> diffusion-results/position-$c-y.tsv
	> diffusion-results/position-$c-z.tsv

	cat results/calcium-ion-diffusion.vec | awk -v r=$((i+0)) '{ if ($1 == r) print $4 }' >> diffusion-results/collisions-$c.tsv

	cat results/calcium-ion-diffusion.vec | awk -v r=$((i+1)) '{ if ($1 == r) print $3"\t"$4 }' >> diffusion-results/collision-position-$c-x.tsv
	cat results/calcium-ion-diffusion.vec | awk -v r=$((i+2)) '{ if ($1 == r) print $3"\t"$4 }' >> diffusion-results/collision-position-$c-y.tsv
	cat results/calcium-ion-diffusion.vec | awk -v r=$((i+3)) '{ if ($1 == r) print $3"\t"$4 }' >> diffusion-results/collision-position-$c-z.tsv

	cat results/calcium-ion-diffusion.vec | awk -v r=$((i+4)) '{ if ($1 == r) print $4 }' >> diffusion-results/position-$c-x.tsv
	cat results/calcium-ion-diffusion.vec | awk -v r=$((i+5)) '{ if ($1 == r) print $4 }' >> diffusion-results/position-$c-y.tsv
	cat results/calcium-ion-diffusion.vec | awk -v r=$((i+6)) '{ if ($1 == r) print $4 }' >> diffusion-results/position-$c-z.tsv
	c=$((c+1))

done