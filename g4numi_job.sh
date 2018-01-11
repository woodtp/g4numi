#!/bin/bash

echo
echo "======== cd to CONDOR_DIR_INPUT ========"
cd $CONDOR_DIR_INPUT

echo
echo "======== ls ========"
ls

echo
echo "======== UNTARRING... ========"
tar xvfz local_install.tar.gz -C ./ > /dev/null

echo
echo "======== Done untarring. ls ========"
ls

echo
echo "======== SETUPS ========"
echo "PROCESS=$PROCESS"
SEED=$((RUN*10000+PROCESS))
echo "SEED=$SEED"
sed -i 's/\${seed}/'$SEED'/g' g4numi.mac

source setup_beamsim.sh
echo "BEAMCONFIG=${BEAMCONFIG}"
echo "PLAYLIST=${PLAYLIST}"
echo "G4NUMIVER=${G4NUMIVER}"

# the outfile is named by the $PROCESS
OUTFILE="g4numi${G4NUMIVER}_${PLAYLIST}_${BEAMCONFIG}_${PROCESS}"
echo $OUTFILE
sed -i 's/\${outfile}/'$OUTFILE'/g' g4numi.mac

echo
echo "======== MACRO CONTENT ========"
cat g4numi.mac

echo
echo "======== EXECUTING g4numi: g4numi g4numi.mac  ========"
g4numi g4numi.mac

echo
echo "Moving output to CONDOR_DIR_G4NUMI"
mv *.root $CONDOR_DIR_G4NUMI
