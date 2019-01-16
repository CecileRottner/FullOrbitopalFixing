#!/bin/bash
rm result.txt

dossier=data/
nom=result_7.txt

n=60
T=48

for sym in 2 ; do
  for id in {1..10}; do
    for met in -1 0 1 5 ; do
      ./mf $met $dossier $n $T 1 3 $sym 1 0 $id $nom
    done
    printf "\n" >> result.txt
  done
done
printf "\n" >> result.txt
printf "\n" >> result.txt

n=30
T=96

for sym in 2 ; do
  for id in {1..10}; do
    for met in -1 0 1 5 ; do
      ./mf $met $dossier $n $T 1 3 $sym 1 0 $id $nom
    done
    printf "\n" >> result.txt
  done
done
printf "\n" >> result.txt
printf "\n" >> result.txt