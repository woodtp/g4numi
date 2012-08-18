#!/usr/bin/env python
import os, re, shutil

#######################################
# This script will rename g4numi
# output to have a fixed-width
# run/subrun number.
# The fixed width lets us select
# a glob of files to use in GENIE.
###########
# example:
#    OLD: /minerva/data/flux/g4numi/v4_test/le100z-200i/g4numiv4_le100z-200i_98_0001.root
#    NEW: /minerva/data/flux/g4numi/v4_test/le100z-200i/g4numiv4_le100z-200i_0098_0001.root
###########
# Brian Tice (tice@physics.rutgers.edu)
# Aug 17, 2012
#######################################

####################################
# Params
####################################
test    = True   #if true don't move, just say what you would do
verbose = True   #print the old and new filenames

g4numiv = "v4"   #what is the g4numiv
fluxdir = "/minerva/data/flux_hadron_samples/flux/g4numi/v4_test"  #what is the top directory for the g4numi files
####################################

for root, dirs, files in os.walk( os.path.join( fluxdir ) ):
  for name in files:
    m = re.search( r'g4numi%s_(\S+)_(\d+)_(\d{4}).root' % g4numiv, name )
    if m:
      beam   = m.group(1)
      run    = int( m.group(2) )
      subrun = int( m.group(3) )
      newname = "g4numi%s_%s_%04d_%04d.root" % (g4numiv, beam, run, subrun )

      oldfile = os.path.join( root, name )
      newfile = os.path.join( root, newname )

      if verbose or test:
        print "============ Moving file ============"
        print "    OLD: %s" % oldfile
        print "    NEW: %s" % newfile

      if test:
        continue

      try:
        shutil.move( oldfile, newfile )
      except Exception, e:
        print "  Problem with move : " + str(e)
