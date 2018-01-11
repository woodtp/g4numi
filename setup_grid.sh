#!/bin/bash

# root
source /cvmfs/minerva.opensciencegrid.org/minerva/beamsim/x86_64/root/bin/thisroot.sh

# misc grid
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup
setup jobsub_client
setup ifdhc

# grid permissions
voms-proxy-destroy # destroy the existing certificates, to get a new with
kx509
voms-proxy-init --rfc --voms=fermilab:/fermilab/minerva/Role=Analysis --noregen
export CPN_LOCK_GROUP=gpcf

# dk2nu
source /cvmfs/nova.opensciencegrid.org/externals/setup
setup dk2nu v01_04_01g -q debug:e9:r5
