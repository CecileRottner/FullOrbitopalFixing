#!/usr/bin/env bash

for file in Results/*500.txt  ; do
	name=${file%500.txt}
	
	name="${name}502.txt"
	mv $file $name
done

