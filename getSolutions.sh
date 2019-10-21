#!/usr/bin/env bash
rm resultat.txt

for n in 20 30 ; do
for T in 48 96; do
  for sym in 2 3; do
    for id in {1..20}; do
	
	for met in 10 11 30 31 40 41 50 51 60 61 70 71; do
	
		cat ResultsRamp/$n_$T_1_3_$sym_0_0_$id_$met.txt >> resultat.txt
	
	printf "\n" >> resultat.txt
	done
  done
done
done