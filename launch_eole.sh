#!/bin/bash
rm Results/*

dossier=data/smaller/

n=20
T=48

  for sym in 2; do
    for id in {1..20}; do
	
	for met in 10 11 30 31 40 41 50 51 60 61 70 71; do
		rm script.sh
		
		echo "#!/usr/bin/env bash" >> "\n" >> script.sh
		echo -n "./mf "  >> script.sh
		echo -n ${met} >> script.sh 
		echo -n " "  >> script.sh 
		echo -n ${dossier} >> script.sh 
		echo -n " "  >> script.sh 
		echo -n ${n} >> script.sh 
		echo -n " "  >> script.sh 
		echo -n ${T} >> script.sh 
		echo -n " "  >> script.sh 
		echo -n 1 >> script.sh 
		echo -n " "  >> script.sh 
		echo -n 3 >> script.sh 
		echo -n " "  >> script.sh 
		echo -n ${sym} >> script.sh 
		echo -n " "  >> script.sh 
		echo -n 0 >> script.sh 
		echo -n " "  >> script.sh 
		echo -n 0 >> script.sh 
		echo -n " "  >> script.sh 
		echo -n ${id} >> script.sh 
   
		chmod +x script.sh
		
        sbatch --exclusive -N 1 --time=12:00:00 --wckey=P11J5:APOGENE script.sh

    done
    printf "\n" >> result.txt
  done
done
printf "\n" >> result.txt
printf "\n" >> result.txt






