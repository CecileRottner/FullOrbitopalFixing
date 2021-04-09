#!/usr/bin/env bash
rm resultat.txt

for n in 60 ; do
for T in 48 96 ; do
  for sym in 2 3 4; do
    for id in {1..20}; do
	#100 101 150 200 300 301 400 401 502 503 600 700
	for met in  100 101 150 151 200 201 300 301 400 401 502 503 600 700; do
	
		fichier=ResultsNoRamp/${n}_${T}_1_3_${sym}_0_0_${id}_${met}.txt
		
		
		    if [ -e ${fichier} ]; then
		
			if [ -s ${fichier}  ]; then
			
				cat ${fichier}  >> resultat.txt
			else
				rm ${fichier} 
				echo "le fichier ${fichier} est vide"
			fi
		    else
			echo "le fichier ${fichier} n'existe pas"
		    fi
	done
	
	printf "\n" >> resultat.txt
	done
  done
done
done
