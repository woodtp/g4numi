#!/bin/bash

#curr=`/bin/pwd`
cd $1
source setup_beamsim.sh
cd $CONDOR_DIR_TUPLES
$1/g4na49 $2


cd $CONDOR_DIR_HISTOS

cp $1/CreateYields.C .
#cp $1/CreateYields_C.so .
files="$CONDOR_DIR_TUPLES/*.root"
echo 'CreateYields.C('$3',"FTFP_BERT","'$files'",'$4')'
root -l -b -n -q 'CreateYields.C('$3',"FTFP_BERT","'$files'",'$4')'
rm ./CreateYields.*


#cp $1/CreateInvXS.C . 
#cp $1/CreateInvXS_C.so . 
#yfile=`ls Yields*.root`
#sfile=`echo $yfile | sed 's/Yields/InvXS/'`
#echo "local area ls"
#ls -l .
#echo 'CreateInvXS.C('$3','$6',"'$yfile'","'$sfile'")'
#root -l -b -n -q 'CreateInvXS.C('$3','$6',"'$yfile'","'$sfile'")'
#rm ./CreateInvXS.*

if [ $5 == 1 ]; then
    echo "keeping the ntuple for transfer back"
elif [ $5 == 0 ]; then
    echo "deleting the output ntuple"
    cd $CONDOR_DIR_TUPLES
    echo `pwd`
    ls ./
    rm ./*.root
fi
