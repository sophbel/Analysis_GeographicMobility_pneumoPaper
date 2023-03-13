#!/bin/bash
variable="$(wc -l ${1} | cut -d' ' -f1)"
newvar="$(expr $variable - 43)"
cat > ${1}_LLs.txt << sed -n "10,${newvar}p" ${1} | cut -d' ' -f2 
