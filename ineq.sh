#!/bin/bash
rm result.txt

dossier=data/Litt_Real/

n=30
T=96

for met in 4; do
  for sym in 4 3 2; do
    for id in {1..20}; do
      ./mf $met $dossier $n $T 1 3 $sym 0 0 $id
    done
    printf "\n" >> result.txt
  done
done
printf "\n" >> result.txt
printf "\n" >> result.txt






