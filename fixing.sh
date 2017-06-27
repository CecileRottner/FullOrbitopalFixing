#!/bin/bash
rm result.txt

dossier=data/Litt_Real/


n=20
T=192

for sym in 4 3; do
  for id in {1..20}; do
      ./mf 0 $dossier $n $T 1 3 $sym 0 0 $id
  done
done
printf "\n" >> result.txt
printf "\n" >> result.txt


n=60
T=48

for sym in 4 3 ; do
  for id in {1..20}; do
      ./mf 0 $dossier $n $T 1 3 $sym 0 0 $id
  done
done
printf "\n" >> result.txt
printf "\n" >> result.txt



n=30
T=96

for sym in 4 3 2 ; do
  for id in {1..20}; do
      ./mf 0 $dossier $n $T 1 3 $sym 0 0 $id
  done
done
printf "\n" >> result.txt
printf "\n" >> result.txt



n=15
T=288

for sym in 3 2 ; do
  for id in {1..20}; do
      ./mf 0 $dossier $n $T 1 3 $sym 0 0 $id
  done
done
