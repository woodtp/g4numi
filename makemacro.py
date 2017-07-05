#!/usr/bin/env python

#$Header: /cvs/projects/numi-beam-sim/numi-beam-sim/g4numi/Attic/makemacro.py,v 1.1.2.14 2017/07/05 16:13:21 bmesserl Exp $

import os, re, sys, getopt,string,optparse
from string import Template

#beamconfig/playlist/targetZpos lookup
beamconfig_dict = { "le010z185i"  : {"minerva1"   : "-44.50",    #tgt2H1 = 9.50
                                     "minerva7"   : "-44.18",    #tgt2H1 = 9.18
                                     "minerva9"   : "-45.40",    #tgt2H1 = 10.40 
                                     "minerva13"  : "-44.17"},   #tgt2H1 = 9.17
                    "le010z-185i" : {"downstream" : "-44.50",    #tgt2H1 = 9.5
                                     "minerva5"   : "-43.85",    #tgt2H1 = 8.85
                                     "minerva10"  : "-44.18"},   #tgt2H1 = 9.18
                    "le010z000i"  : {"minerva6"   : "-44.18"},   #tgt2H1 = 9.18
                    "le100z200i"  : {"minerva2"   : "-134.57",   #tgt2H1 = 99.57
                                     "minerva11"  : "-134.17"},  #tgt2H1 = 99.17
                    "le100z-200i" : {"minerva3"   : "-134.57",   #tgt2H1 = 99.57
                                     "minerva12"  : "-134.17"},  #tgt2H1 = 99.17
                    "le250z200i"  : {"minerva4"   : "-284.57",   #tgt2H1 = 249.57
                                     "minerva8"   : "-285.09"}   #tgt2H1 = 250.09
                  }

#Defaults
HORN1_POSITION_X  = 0      #cm
HORN1_POSITION_Y  = 0      #cm

HORN2_POSITION_X  = 0      #cm

BEAM_POSITION_X   = 0      #m
BEAM_POSITION_Y   = 0      #m

BEAM_SPOTSIZE_X   = 1.4    #mm (ME!)
BEAM_SPOTSIZE_Y   = 1.4    #mm (ME!)

TARGET_POSITION_X = 0.0    #cm
TARGET_POSITION_Y = 0.0    #cm
TARGET_POSITION_Z = -143.3 #cm (ME!)

HORN_WATER_MM    = 1       #mm 
POT              = 400000
RUN              = 1
TEMPLATE         = "{0}/macros/template_ME.mac".format(os.getenv("BEAMSIM"))

TARGET_WATER_CM = 3      #cm
BEAMCONFIG        = "me000z200i"
PLAYLIST          = "minervame"

def get_options():
  parser        = optparse.OptionParser(usage="usage: %prog [options]")
  horn_group    = optparse.OptionGroup(parser, "Horn Options")
  target_group  = optparse.OptionGroup(parser, "Target Options")
  beam_group    = optparse.OptionGroup(parser, "Beam Options")
  job_group     = optparse.OptionGroup(parser, "Job Options")
  old_group     = optparse.OptionGroup(parser, "Old Options")

  horn_group.add_option('--do_horn1_old_geometry',   action="store_true", 
        help="The 'old' horn1 geom (formerly known as 'alternate') is now default. " 
             "Use this option to use the old geometry.")
  horn_group.add_option('--do_horn1_fine_segmentation', action="store_true", 
        help="Works for old and new horn1.")

  horn_group.add_option('--horn1_position_X', 
        help="horn 1 transverse offset (_X0 position). In cm.")
  horn_group.add_option('--horn1_position_Y', 
        help="horn 1 vertical offset (_Y0 position). In cm.")
  horn_group.add_option('--horn2_position_X', 
        help="horn 2 transverse offset (_X0 position). In cm.")

  horn_group.add_option('--horn_water_mm', 
        help="mm water layer on horn. ")

  target_group.add_option('--target_position_X', 
          help="target horizontal position.")
  target_group.add_option('--target_position_Y', 
          help="target vertical position.")
  target_group.add_option('--target_position_Z', 
          help="target longitudinal position."
               "NOTE: this will overwrite target Z set by beamconfig.")

  beam_group.add_option('--beam_position_X', 
        help="beam horizontal position.")
  beam_group.add_option('--beam_position_Y', 
        help="beam vertical position.")
  beam_group.add_option('--beam_spotsize_X', 
        help="beam horizontal spot size.")
  beam_group.add_option('--beam_spotsize_Y', 
        help="beam vertical spot size.")

  job_group.add_option('--pot',
        help="default number of protons on target to simulate")
  job_group.add_option('--run',
        help="Run number.")
  job_group.add_option('--template',
        help='specify template macro')
  job_group.add_option('--filetag',
        help="doesn't come with the underscore")

  old_group.add_option('--no_importance_weighting', action="store_true", 
        help="No importance weighting." )
  old_group.add_option('--seed',
        help="Random seed used by g4numi.")
  old_group.add_option('--target_water_cm',
        help="simulate water in the target(cm)." )
  old_group.add_option('--playlist',
        help='correct the Z target posision per playlist (LE)')
  old_group.add_option('--beamconfig', 
        help="Configure the neutrino beam. "
              "This sets the TargetZ0 and Baffle Z0 and HornCurrent. "
              "Must be in the form 'LE#z#i' or 'le#z#i', where # is any number. "
              "Examples are le010z185i, LE025.3z-200i, LE250z185.6i....etc.")

  parser.add_option_group(horn_group)
  parser.add_option_group(beam_group)
  parser.add_option_group(target_group)
  parser.add_option_group(job_group)
  parser.add_option_group(old_group)

  options, remainder = parser.parse_args()
  return options

def main():
  options = get_options()

  #set defaults

  #use beamconfig and playlist to make a first determination of target position and horn current
  #at the moment, this is the only way to change horn current
  #but the dedicated target z position option will overwrite this beamconfig option
  beamconfig                 = BEAMCONFIG        if not options.beamconfig                 else options.beamconfig
  playlist                   = PLAYLIST.lower()  if not options.playlist                   else options.playlist.lower()

  # different beamconfig has been specified
  if beamconfig != BEAMCONFIG:
    try:
      target_position_Z  = beamconfig_dict[beamconfig][playlist] if not options.target_position_Z else options.target_position_Z
    except:
      print sys.exit("Error! beamconfig-playlist pair not found!")
  else:
    target_position_Z        = TARGET_POSITION_Z if not options.target_position_Z          else options.target_position_Z #cm

  #horn
  do_horn1_fine_segmentation = False             if not options.do_horn1_fine_segmentation else True
  do_horn1_new_geometry      = True              if not options.do_horn1_old_geometry      else False
  horn_water_mm              = HORN_WATER_MM     if not options.horn_water_mm              else options.horn_water_mm #mm
  horn1_position_X           = HORN1_POSITION_X  if not options.horn1_position_X           else options.horn1_position_X #cm
  horn1_position_Y           = HORN1_POSITION_Y  if not options.horn1_position_Y           else options.horn1_position_Y #cm
  horn2_position_X           = HORN2_POSITION_X  if not options.horn2_position_X           else options.horn2_position_X #cm

  #beam size/position
  beam_position_X            = BEAM_POSITION_X   if not options.beam_position_X            else options.beam_position_X #m
  beam_position_Y            = BEAM_POSITION_Y   if not options.beam_position_Y            else options.beam_position_Y #m
  beam_spotsize_X            = BEAM_SPOTSIZE_X   if not options.beam_spotsize_X            else options.beam_spotsize_X #mm
  beam_spotsize_Y            = BEAM_SPOTSIZE_Y   if not options.beam_spotsize_Y            else options.beam_spotsize_Y #mm

  #target
  target_position_X          = TARGET_POSITION_X if not options.target_position_X          else options.target_position_X #cm
  target_position_Y          = TARGET_POSITION_Y if not options.target_position_Y          else options.target_position_Y #cm
  do_target_water            = False             if not options.target_water_cm            else True
  target_water_cm            = TARGET_WATER_CM   if not options.target_water_cm            else options.target_water_cm 

  #job
  filetag                    = ""                if not options.filetag                    else '_'+options.filetag
  pot                        = POT               if not options.pot                        else options.pot
  run                        = RUN               if not options.run                        else options.run
  template                   = TEMPLATE          if not options.template                   else options.template

  #misc
  do_importance_weighting    = True              if not options.no_importance_weighting    else False
  seed                       = run               if not options.seed                       else options.seed 

  outfile = "g4numi{0}_{1}_{2}".format(os.getenv("G4NUMIVER"), playlist, beamconfig)
  outfile = outfile + filetag

  filestring=open(template,'r').read()
  t=string.Template(filestring)
  print t.substitute({'do_horn1_new_geometry': do_horn1_new_geometry,
                      'do_horn1_fine_segmentation': do_horn1_fine_segmentation,
                      'horn1_position_X': horn1_position_X,
                      'horn1_position_Y': horn1_position_Y,
                      'horn2_position_X': horn2_position_X,
                      'horn_water_mm': horn_water_mm,
                      'beam_position_X': beam_position_X,
                      'beam_position_Y': beam_position_Y,
                      'beam_spotsize_X': beam_spotsize_X,
                      'beam_spotsize_Y': beam_spotsize_Y,
                      'target_position_X': target_position_X,
                      'target_position_Y': target_position_Y,
                      'target_position_Z': target_position_Z,
                      'target_water_cm': target_water_cm,
                      'do_target_water' : do_target_water,
                      'pot':pot,
                      'run':run,
                      'beamconfig':beamconfig,
                      'playlist': playlist,
                      'do_importance_weighting': do_importance_weighting,
                      'outfile':outfile,
                      'seed':seed})

if __name__ == "__main__":
  sys.exit(main())
