#!/usr/bin/env bash
rm resultat.txt

for n in 20 60 ; do
for T in 48 96; do
  for sym in 2 3 4; do
    for id in {1..20}; do
	
	for met in 52 53; do
	
		fichier=ResultsRamp/${n}_${T}_1_3_${sym}_0_0_${id}_${met}.txt
		
		
		    if [ -e ${fichier} ]; then
		
			if [ -s ${fichier}  ]; then
			
				head -1 ${fichier}  >> resultat.txt
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
