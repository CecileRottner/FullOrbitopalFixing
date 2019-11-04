#!/usr/bin/env bash
rm resultat.txt

for n in 20 30 ; do
for T in 48 96; do
  for sym in 2 3; do
    for id in {1..20}; do
	
	for met in 10 11 20 21 30 31 40 41 50 51 60 70; do
	
		fichier=ResultsRamp/${n}_${T}_1_3_${sym}_0_0_${id}_${met}.txt
		
		
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
