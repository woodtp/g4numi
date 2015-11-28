#!/bin/bash
echo $X509_USER_PROXY
# must use -d option to minerva_jobsub 
# e.g.,  -d G4NUMI /minerva/data/flux/g4numi/version/beamconfig/
echo "CONDOR_DIR_G4NUMI = ${CONDOR_DIR_G4NUMI}"
cd $CONDOR_DIR_G4NUMI

#setup g4numi - sets BEAMSIM,
# /minerva/app/users/$USER/cmtuser/NumiAna/numisoft/g4numi
# must pass in the G4NUMIAREA via -e G4NUMIAREA
echo "G4NUMIAREA = $G4NUMIAREA"
WORKDIR=$G4NUMIAREA
source $WORKDIR/setup_beamsim.sh

setup_beamsim
#create the macro file

echo "BEAMCONFIG=$BEAMCONFIG"
echo "DOWATER=$DOWATER"
echo "DONEWHORN1=$DONEWHORN1"
echo "WATERCM=$WATERCM"
echo "DOIMPWT=$DOIMPWT"
echo "POT=$POT"
echo "RUN=$RUN"
echo "PROCESS=$PROCESS"
SEED=$((RUN*10000+PROCESS))
echo "SEED=$SEED"

#check for ME input:
metag1="me"
metag2="ME"
templatename="template.mac"
case $BEAMCONFIG in 
    *"$metag1"*)templatename="template_ME.mac";;
    *"$metag2"*)templatename="template_ME.mac";;
esac

#check for LE target position input (all cases that it is applicable):

case $BEAMCONFIG in 
    "le010z185i"|"LE010z185i")TGTPOS="-z $PLAYLIST";;
    "le010z-185i"|"LE010z-185i")TGTPOS="-z $PLAYLIST";;
    "le010z000i"|"LE010z000i")TGTPOS="$-z $PLAYLIST";;
    "le100z200i"|"LE100z200i")TGTPOS="$-z $PLAYLIST";;
    "le100z-200i"|"LE100z-200i")TGTPOS="-z $PLAYLIST";;
    "le250z200i"|"LE250z200i")TGTPOS="-z $PLAYLIST";;
esac
echo "TGTPOS=$TGTPOS"

IFLAG=""
if [ "$DOIMPWT" = "false" ]; then
    IFLAG="-i"
fi

MFLAG=""
if [ "$DONEWHORN1" = "true" ]; then
    MFLAG="-m"
fi

if [ "$DOWATER" = "false" ]; then
    
    $WORKDIR/makemacro.py -t $WORKDIR/macros/$templatename -b $BEAMCONFIG ${MFLAG} -p $POT -r $RUN -n $PROCESS -s $SEED -z $PLAYLIST $IFLAG > g4numi.mac

else 

    $WORKDIR/makemacro.py -t $WORKDIR/macros/$templatename -b $BEAMCONFIG -w ${MFLAG} -L $WATERCM -p $POT -r $RUN -n $PROCESS -s $SEED -z $PLAYLIST $IFLAG > g4numi.mac
	
fi

cat g4numi.mac
$WORKDIR/g4numi g4numi.mac
rm g4numi.mac
