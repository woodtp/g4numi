#!/bin/bash

moms="158 120 110 100 90 80 70 60 50 40 31 20 12"
runlow="1"
runhigh="500"
for m in $moms; do
./process_na49Grid.pl ${m} ${runlow} ${runhigh} -p FTFP_BERT -o /minerva/data/users/kordosky/hp
done
