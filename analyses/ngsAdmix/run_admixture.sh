#!/bin/bash

#parallel NGSadmixConv.sh beagle_test.beagle.gz {1} ::: {3..9}

NGSCONV="NGSadmixConv.sh"

bgl=muskox_muskoxmap_whitelist.beagle.gz

seq 5 13 | xargs -n1 -P2 $NGSCONV $bgl

#for k in `seq 2 6`
#do
#        $NGSCONV $bgl $k
#done
