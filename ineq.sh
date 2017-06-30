#!/bin/bash
rm result.txt

dossier=data/FixedSym/


met=0
n=60
T=48

for sym in 4 3 2 ; do
  for id in {1..20}; do
      ./mf $met $dossier $n $T 1 3 $sym 0 0 $id
  done
done
printf "\n" >> result.txt
printf "\n" >> result.txt





