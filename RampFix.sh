#!/bin/bash
rm result.txt

dossier=data/newdata/

nom=result_7.txt

n=60
T=96

for sym in 2 ; do
  for id in {1..20}; do
    for met in -1 0 1 4 44 ; do
      ./mf $met $dossier $n $T 1 3 $sym 0 0 $id $nom
    done
    printf "\n" >> result.txt
  done
done
printf "\n" >> result.txt
printf "\n" >> result.txt
