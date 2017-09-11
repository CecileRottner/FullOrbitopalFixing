#!/bin/bash
rm result.txt

dossier=data/Litt_Real/

n=20
T=192

for met in 2; do
  for sym in 2; do
    for id in {2..5}; do
      ./mf $met $dossier $n $T 1 3 $sym 0 0 $id
    done
  done
done
printf "\n" >> result.txt
printf "\n" >> result.txt






