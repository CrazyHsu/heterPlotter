#!/bin/bash
<<BLOCK
File name: test.sh
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2018-10-26 18:37:48
Last modified: 2019-08-13 00:51:19
BLOCK

python heterPlotter.py -g data/chr3.gtf -i 3:158397000-158598000 -hp data/test.hapmap -gw data/test.gwas -a data/anno.txt -o test.png
rm *.idx *.num *.pect *.tmp
