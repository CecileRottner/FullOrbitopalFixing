#!/bin/bash


dossier=data/newdata/

#sym=2, id=6, met=10
#sym=3, id= 2 13 15 16 17 18 19 20, met=10 11 50 51

for n in 60 ; do
for T in 96 ; do
  for sym in 2 3 4; do
    for id in {1..20}; do
	
	for met in 100 101 150 151 200 201 300 301 400 401 500 501 600 700 ; do
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
		
        sbatch --exclusive -N 1 --time=01:20:00 --wckey=P11J5:APOGENE script.sh

    done
  done
done
done
done
