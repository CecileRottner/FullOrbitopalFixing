#!/bin/bash

dossier=data/newdata/

nom=result_2.txt



n=60
T=96

for sym in 4 3 2 ; do
  for id in {1..20}; do
    for met in 2 ; do
      ./mf $met $dossier $n $T 1 3 $sym 0 0 $id $nom
    done
    printf "\n" >> $nom
  done
done
printf "\n" >> $nom
printf "\n" >> $nom
