source ../environment.sh
DIR='../data/'$DATA'/'
EDGE=$DIR$DATA'.edge'
COST=$DIR$DATA'.cost'
UI=$DIR$DATA'.user.interest'
KT=$DIR$DATA'.keyword'
IL=$DIR$DATA'.il'
QUERY=$DIR$DATA'.query'


# USET=$DIR$DATA'.nb'
# bound1=1000
# echo 'offline'
# make clean && make offline
# echo "$EDGE $COST $bound1 $USET"
# ./offline -nb $EDGE $COST $bound1 $USET

USET=$DIR$DATA'.knb'
bound1=1000
echo 'offline'
make clean && make offline
echo "$EDGE $COST $IL $KT $bound1 $USET"
./offline -knb $EDGE $COST $KT $UI $IL $bound1 $USET

