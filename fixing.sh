#!/bin/bash
rm result.txt

dossier=data/newdata/

nom=result_1.txt

n=60
T=48

for sym in 5 ; do
  for id in {1..20}; do
      ./mf 1 $dossier $n $T 1 3 $sym 0 0 $id $nom
  done
done
printf "\n" >> result.txt
printf "\n" >> result.txt


n=80
T=48

for sym in 6 5 4 ; do
  for id in {1..20}; do
      ./mf 1 $dossier $n $T 1 3 $sym 0 0 $id
  done
done
printf "\n" >> result.txt
printf "\n" >> result.txt


n=60
T=96

for sym in 4 3 2 ; do
  for id in {1..20}; do
      ./mf 1 $dossier $n $T 1 3 $sym 0 0 $id
  done
done
printf "\n" >> result.txt
printf "\n" >> result.txt


