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
#    OLD: /minerva/data/flux/g4numi/v5/le100z-200i/g4numiv5_minerva1_le010z185i_0111_0005.root
#    NEW: /minerva/data/flux/g4numi/v5/le100z-200i/g4numiv5_dk2nu_minerva1_le010z185i_0111_0005.root
###########
# Brian Tice (tice@physics.rutgers.edu)
# Aug 17, 2012
#######################################
# This script was slightly modifed to add dk2nu new format name
#(Leo Aliaga, Feb 2015)
#######################################



test    = True
verbose = True  

g4numiv = "v5" 
fluxdir = "/minerva/data/users/laliaga/flux/le100z200i"

####################################

for root, dirs, files in os.walk( os.path.join( fluxdir ) ):
  for name in files:
    m = re.search( r'g4numi%s_(\S+)_(\S+)_(\d+)_(\d{4}).root' % g4numiv, name )
    if m:
      pl     = m.group(1)
      beam   = m.group(2)
      run    = int( m.group(3) )
      subrun = int( m.group(4) )
      newname = "g4numi%s_dk2nu_%s_%s_%04d_%04d.root" % (g4numiv, pl , beam, run, subrun )

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

        
