source ../environment.sh

DIR='../data/'$DATA'/'
EDGE=$DIR$DATA'.edge'
COST=$DIR$DATA'.cost'
UI=$DIR$DATA'.user.interest'
KT=$DIR$DATA'.keyword'
NB=$DIR$DATA'.knb'
IL=$DIR$DATA'.il'
QUERY=$DIR$DATA'.query'

bound1=1000
bound2=10

make clean && make main
mv main exp_best

BESTFIRST='../figs/'$DATA'/keywordBudget'
rm  $BESTFIRST
touch $BESTFIRST
for B in 20
do
    echo "./experiments -b $EDGES $COST $UI $KT $IL $NB $QUERY"
    ./exp_best -b $EDGE $COST $UI $KT $IL $NB $QUERY $B $bound1 $bound2 >>$BESTFIRST
done

make clean
rm exp_best



