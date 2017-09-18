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
export HOME=${TMPDIR}

#create the macro file

echo "BEAMCONFIG=$BEAMCONFIG"
echo "DO_TARGET_WATER=$DO_TARGET_WATER"
echo "DO_HORN1_OLD_GEOMETRY=$DO_HORN1_OLD_GEOMETRY"
echo "DO_HORN1_FINE_SEGMENTATION=$DO_HORN1_FINE_SEGMENTATION"
echo "IMPORTANCE_WEIGHTING=$IMPORTANCE_WEIGHTING"
echo "TARGET_WATER_CM=$TARGET_WATER_CM"
echo "POT=$POT"
echo "RUN=$RUN"
echo "PROCESS=$PROCESS"
SEED=$((RUN*10000+PROCESS))
echo "SEED=$SEED"
echo "HORN1_POSITION_X=$HORN1_POSITION_X"
echo "HORN1_POSITION_Y=$HORN1_POSITION_Y"
echo "HORN1_POSITION_Z=$HORN1_POSITION_Z"
echo "HORN2_POSITION_X=$HORN2_POSITION_X"
echo "HORN2_POSITION_Y=$HORN2_POSITION_Y"
echo "TARGET_POSITION_X=$TARGET_POSITION_X"
echo "TARGET_POSITION_Y=$TARGET_POSITION_Y"
echo "TARGET_POSITION_Z=$TARGET_POSITION_Z"
echo "HORN_WATER_MM=$HORN_WATER_MM"
echo "BEAM_POSITION_X=$BEAM_POSITION_X"
echo "BEAM_POSITION_Y=$BEAM_POSITION_Y"
echo "BEAM_SPOTSIZE_X=$BEAM_SPOTSIZE_X"
echo "BEAM_SPOTSIZE_Y=$BEAM_SPOTSIZE_Y"

#check for ME input:
metag1="me"
metag2="ME"
templatename="template.mac"
case $BEAMCONFIG in 
    *"$metag1"*)templatename="template_ME.mac";;
    *"$metag2"*)templatename="template_ME.mac";;
esac

echo "templatename = $templatename"

#check for LE target position input (all cases that it is applicable):

#beamconfig setting discontinued in ME
#target position is set manually in g4numi_grid_submit_ME.sh
case $BEAMCONFIG in 
    "le010z185i"|"LE010z185i")TGTPOS="--playlist $PLAYLIST";;
    "le010z-185i"|"LE010z-185i")TGTPOS="--playlist $PLAYLIST";;
    "le010z000i"|"LE010z000i")TGTPOS="$--playlist $PLAYLIST";;
    "le100z200i"|"LE100z200i")TGTPOS="$--playlist $PLAYLIST";;
    "le100z-200i"|"LE100z-200i")TGTPOS="--playlist $PLAYLIST";;
    "le250z200i"|"LE250z200i")TGTPOS="--playlist $PLAYLIST";;
esac
echo "TGTPOS=$TGTPOS"

# set some bool flags for the make macro py script
DO_HORN1_OLD_GEOMETRY_FLAG=""
if [ "$DO_HORN1_OLD_GEOMETRY" = "true" ]; then
    DO_HORN1_OLD_GEOMETRY_FLAG="--do_horn1_old_geometry"
fi

DO_HORN1_FINE_SEGMENTATION_FLAG=""
if [ "$DO_HORN1_FINE_SEGMENTATION" = "true" ]; then
    DO_HORN1_FINE_SEGMENTATION_FLAG="--do_horn1_fine_segmentation"
fi

NO_IMPORTANCE_WEIGHTING_FLAG=""
if [ "$IMPORTANCE_WEIGHTING" = "false" ]; then
    NO_IMPORTANCE_WEIGHTING_FLAG="--no_importance_weighting "
fi

TARGET_WATER_FLAG=""
if [ "$DO_TARGET_WATER" = "true" ]; then
    TARGET_WATER_FLAG="--target_water_cm=$TARGET_WATER_CM"
fi

# make the macro
$WORKDIR/makemacro.py \
  $DO_HORN1_OLD_GEOMETRY_FLAG \
  $DO_HORN1_FINE_SEGMENTATION_FLAG \
  --horn1_position_X $HORN1_POSITION_X \
  --horn1_position_Y $HORN1_POSITION_Y \
  --horn1_position_Z $HORN1_POSITION_Z \
  --horn2_position_X $HORN2_POSITION_X \
  --horn2_position_Y $HORN2_POSITION_Y \
  --horn_water_mm $HORN_WATER_MM \
  --target_position_X $TARGET_POSITION_X \
  --target_position_Y $TARGET_POSITION_Y \
  --target_position_Z $TARGET_POSITION_Z \
  --beam_position_X $BEAM_POSITION_X \
  --beam_position_Y $BEAM_POSITION_Y \
  --beam_spotsize_X $BEAM_SPOTSIZE_X \
  --beam_spotsize_Y $BEAM_SPOTSIZE_Y \
  --pot $POT \
  --run $RUN \
  --template $WORKDIR/macros/$templatename \
  --filetag $PROCESS \
  $NO_IMPORTANCE_WEIGHTING_FLAG \
  --seed $SEED \
  $TARGET_WATER_FLAG \
  --playlist $PLAYLIST \
  --beamconfig $BEAMCONFIG \
  > g4numi.mac

cat g4numi.mac
$WORKDIR/g4numi g4numi.mac
rm g4numi.mac
