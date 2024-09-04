#!/bin/bash

echo
echo "======== cd to CONDOR_DIR_INPUT ========"
cd $CONDOR_DIR_INPUT

echo
echo "======== ls ========"
ls

echo
echo "======== bootstrap a number of UPS areas  ========"
export   GPRDCVMFS=/cvmfs/fermilab.opensciencegrid.org/products/genie
source ${GPRDCVMFS}/bootstrap_genie_ups.sh

echo
echo "======== UNTARRING... ========"
# tar xvfz local_products.tar.gz -C ./ > /dev/null

echo
echo "======== Done untarring. ls ========"
ls -l ${INPUT_TAR_DIR_LOCAL}/local_products

echo
echo "======== setup our version of g4num as a UPS product  ========"
export PRODUCTS=${PRODUCTS}:${INPUT_TAR_DIR_LOCAL}/local_products
ups list -aK+ g4numi
xyzzy=`ups list -aK+ g4numi | tr -d '"' `
version=`echo $xyzzy | cut -d' ' -f2`
qualifier=`echo $xyzzy | cut -d' ' -f4`

setup g4numi $version -q $qualifier




echo
echo "======== UPDATE MACRO WITH RUN NUMBER ========"
SEED=$(date +%Y%m%d%H%M%S)
sed -i 's/\${seed}/'$SEED'/g' g4numi.mac
OUTFILE="g4numi${G4NUMIVER}_${PLAYLIST}_${BEAMCONFIG}_${PROCESS}_${SEED}"
sed -i 's/\${outfile}/'$OUTFILE'/g' g4numi.mac

echo "BEAMCONFIG=${BEAMCONFIG}"
echo "PLAYLIST=${PLAYLIST}"
echo "G4NUMIVER=${G4NUMIVER}"
echo "PROCESS=$PROCESS"
echo "SEED=$SEED"
echo "OUTFILE=$OUTFILE"

echo
echo "======== MACRO CONTENT ========"
cat g4numi.mac

echo
echo "======== EXECUTING g4numi ========"
echo "g4numi g4numi.mac"
g4numi g4numi.mac

echo
#echo "Moving output to CONDOR_DIR_G4NUMI"
#echo "CONDOR_DIR_G4NUMI=$CONDOR_DIR_G4NUMI"
#mv *.root $CONDOR_DIR_G4NUMI

echo "Using ifdh cp to copy output to OUTDIR"
setup ifdhc
for f in *.root ; do
  echo ifdh cp $f ${OUTDIR}/$f
       ifdh cp $f ${OUTDIR}/$f
  echo status=$?
done


# end-of-script
