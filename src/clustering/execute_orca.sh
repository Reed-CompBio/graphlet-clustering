#!/bin/bash

type=$1
upto=$2
path=${3?Error: No path given}

orca_path=$(pwd)/orca/

cd $path

for i in *.orca
do	
	echo "running: ${orca_path}'orca.exe' ${type} ${upto} ${i} ${i}'.ocount'";
	${orca_path}'orca.exe' ${type} ${upto} ${i} ${i}'.ocount'
done
