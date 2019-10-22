#!/bin/bash
rm Results/*

dossier=data/smaller/


for n in 20  ; do
for T in 48 ; do
  for sym in 2; do
    for id in {1..1}; do
	
	for met in 60 ; do
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
  done
done
done
done


