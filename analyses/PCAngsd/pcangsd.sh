#!/bin/bash

beagle="muskox_muskoxmap_whitelist.beagle.gz"
samples="samples.txt"

python3.5 pcangsd.py -beagle $beagle -sites_save -e 4 -tree -tree_samples $samples \
-o pcangsd -threads 12 > pcangsd.log


#### with Sib
beagle="muskox_muskoxmap_whitelist_wSib.beagle.gz"
samples="samples_wSib.txt"

python3.5 pcangsd.py -beagle $beagle -sites_save -e 4 -tree -tree_samples $samples \
-o pcangsd_wSib -threads 12 > pcangsd_wSib.log


#### sheepmap
beagle="muskox_sheepmap_whitelist.beagle.gz"
samples="samples.txt"

python3.5 pcangsd.py -beagle $beagle -sites_save -e 4 -tree -tree_samples $samples \
-o pcangsd_sheepmap -threads 12 > pcangsd_sheepmap.log

#### sheepmap w Sib
beagle="muskox_sheepmap_whitelist_wSib.beagle.gz"
samples="samples_wSib.txt"

python3.5 pcangsd.py -beagle $beagle -sites_save -e 4 -tree -tree_samples $samples \
-o pcangsd_sheepmap_wSib -threads 12 > pcangsd_sheepmap_wSib.log
