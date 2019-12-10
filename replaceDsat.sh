#!/usr/bin/env bash

for file in Results/*500.txt  ; do
	name=${file%500.txt}
	
	name="${name}502.txt"
	mv $file $name
done

for file in Results/*501.txt  ; do
	name=${file%501.txt}
	
	name="${name}503.txt"
	mv $file $name
done

for file in Results/*50.txt  ; do
	name=${file%50.txt}
	
	name="${name}54.txt"
	mv $file $name
done

for file in Results/*51.txt  ; do
	name=${file%51.txt}
	
	name="${name}55.txt"
	mv $file $name
done

for file in ResultsRamp/60_48*52.txt  ; do
	name=${file%52.txt}
	
	name="${name}54.txt"
	mv $file $name
done

for file in ResultsRamp/60_48*53.txt  ; do
	name=${file%53.txt}
	
	name="${name}55.txt"
	mv $file $name
done


