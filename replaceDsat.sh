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
	
	name="${name}52.txt"
	mv $file $name
done

for file in Results/*51.txt  ; do
	name=${file%51.txt}
	
	name="${name}53.txt"
	mv $file $name
done

