#!/bin/bash
rm result.txt

dossier=data/smaller/

n=20
T=48

  for sym in 2; do
    for id in {1..20}; do
	
	for met in 20 21; do
      ./mf $met $dossier $n $T 1 3 $sym 0 0 $id
    done
    printf "\n" >> result.txt
  done
done
printf "\n" >> result.txt
printf "\n" >> result.txt






