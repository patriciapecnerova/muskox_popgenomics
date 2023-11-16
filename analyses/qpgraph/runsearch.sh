

RSCRIPT="/home/genis/software/R-4.0.2/bin/Rscript"
SEARCHER="/home/genis/github/QCSeq/std_analyses/qpgraph/scriptsfuns/graphExplorer.R"

usepops="Outgroup,SIB,CaME,CaMW,CaIS,CaIN,GrN"

f2dir="/home/genis/other/muskox/qpgraph2_sheepmap/f2_statistics"
nruns=100
outpop="Outgroup"


for nadmix in 0 1 2 3
do
    outpre=/home/genis/other/muskox/qpgraph2_sheepmap/automatic_search/results/muskox_graph_n${nadmix}
    echo "$RSCRIPT $SEARCHER --f2dir $f2dir --outprefix $outpre --outpop $outpop --usepops $usepops --nruns $nruns --nadmix $nadmix"
done | xargs -L1 -I CMD -P20 nice bash -c CMD
