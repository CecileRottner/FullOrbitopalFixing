#!/bin/bash
rm result.txt

dossier=data/smaller/

n=20
T=48

  for sym in 2; do
    for id in {1..20}; do
	
	for met in 10 11 30 31 40 41 50 51 60 61 70 71; do
      ./mf $met $dossier $n $T 1 3 $sym 0 0 $id
    done
    printf "\n" >> result.txt
  done
done
printf "\n" >> result.txt
printf "\n" >> result.txt






