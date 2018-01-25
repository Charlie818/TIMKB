source ../environment.sh

DIR='../data/'$DATA'/'
EDGE=$DIR$DATA'.edge'
COST=$DIR$DATA'.cost'
UI=$DIR$DATA'.user.interest'
KT=$DIR$DATA'.keyword'
NB=$DIR$DATA'.nb'
IL=$DIR$DATA'.il'
QUERY=$DIR$DATA'.query'

bound1=1000
bound2=10

make clean && make main
mv main exp_neigh

BESTFIRST='../figs/'$DATA'/neighbor'
rm  $BESTFIRST
touch $BESTFIRST
for B in 20
do
    echo "./experiments -n $EDGES $COST $UI $KT $IL $NB $QUERY"
    ./exp_neigh -n $EDGE $COST $UI $KT $IL $NB $QUERY $B $bound1 $bound2 >>$BESTFIRST
done
make clean
rm exp_neigh



