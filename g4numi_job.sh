#!/bin/bash
echo $X509_USER_PROXY
# must use -d option to minerva_jobsub 
# e.g.,  -d G4NUMI /minerva/data/flux/g4numi/version/beamconfig/
#echo "CONDOR_DIR_G4NUMI" $CONDOR_DIR_G4NUMI
cd $CONDOR_DIR_G4NUMI

#setup g4numi - sets BEAMSIM,
 
WORKDIR=/minerva/app/users/laliaga/cmtuser/Minerva_v10r6/NumiAna/numisoft/g4numi
source $WORKDIR/setup_beamsim.sh

setup_beamsim -w 1
#create the macro file

echo "BEAMCONFIG=$BEAMCONFIG"
echo "DOWATER=$DOWATER"
echo "WATERCM=$WATERCM"
echo "POT=$POT"
echo "RUN=$RUN"
echo "PROCESS=$PROCESS"
SEED=$((RUN*10000+PROCESS))
echo "SEED=$SEED"
if [ "$DOWATER" = "false" ]; then
    
    $WORKDIR/makemacro.py -t $WORKDIR/macros/template.mac -b $BEAMCONFIG -p $POT -r $RUN -n $PROCESS -s $SEED > g4numi.mac

else 

    $WORKDIR/makemacro.py -t $WORKDIR/macros/template.mac -b $BEAMCONFIG -w -L $WATERCM -p $POT -r $RUN -n $PROCESS -s $SEED > g4numi.mac
	
fi

cat g4numi.mac
$WORKDIR/g4numi g4numi.mac
rm g4numi.mac
