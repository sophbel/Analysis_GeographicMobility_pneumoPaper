#!/bin/bash
variable="$(wc -l ${1} | cut -d' ' -f1)"
newvar="$(expr $variable - 43)"
sed -n "10,${newvar}p" ${1} | cut -d' ' -f2 > ${1}_LLs
