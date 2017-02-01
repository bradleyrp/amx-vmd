{

'vmd_protein':{
#####
####
###
##
#
'quick':"view_routine()",
'extensions':[],'tags':[],'imports':['@vmd'],'params':None,
'settings':"""

step: v01-look
video name: video
view mode: video

#---video settings
max frames: 300         # aim for this many frames in the unbroken trajectory
video size: 30          # size in megabytes (requires ffmpeg 2-pass encoding)
duration: 0.0           # desired duration (requires ffmpeg)

#---resolution
viewbox: (800,800)      # pixels in the VMD viewer (must match the resolution proportions)
resolution: (2400,1800) # snapshot and video resolution
which view: xview       # camera direction
scale: scale by 1.75    # zoom factor (may require some tuning)

#---recipes set the aesthetics of the rendering
recipe collection: video aamd atomistic bilayer protein

#---align the protein
backbone align name: "name CA"

#---custom selections (see lib_vmdmake for details)
selections:| [{'basic_residues':'protein and (resname HIS or resname ARG or resname LYS)',
  'smooth':True,'style':'licorice','goodsell':True},
  ][1:0]
#---^^^ customizations turned off if the list is sliced over 1:0"!
#---note: do not use smooth on some items but not others (induces a non-physical asynchronicity)
#---note: SEE lib_vmdmake.py and the vmdmake documentation for more details

"""},

'vmd_cgmd_bilayer':{
#####
####
###
##
#
'quick':'view_routine()',
'extensions':[],'tags':[],'imports':['@vmd'],'params':None,
'settings':"""

step: v01-look
video name: video
view mode: video

#---video settings
max frames: 300         # aim for this many frames in the unbroken trajectory
video size: 30          # size in megabytes (requires ffmpeg 2-pass encoding)
duration: 0.0           # desired duration (requires ffmpeg)

#---resolution
viewbox: (800,800)      # pixels in the VMD viewer (must match the resolution proportions)
resolution: (2400,1800) # snapshot and video resolution
which view: xview       # camera direction
scale: scale by 1.75    # zoom factor (may require some tuning)

#---recipes set the aesthetics of the rendering
recipe collection: video cgmd bilayer protein

#---custom selections (see lib_vmdmake for details)
selections:| [{'basic_residues':'protein and (resname HIS or resname ARG or resname LYS)',
  'smooth':True,'style':'licorice','goodsell':True},
  ][1:0]
#---^^^ customizations turned off if the list is sliced over 1:0"!
#---note: do not use smooth on some items but not others (induces a non-physical asynchronicity)
#---note: SEE lib_vmdmake.py and the vmdmake documentation for more details

#---COARSE-GRAIN VISUALS
bonder: true
cg_bonds: @cg_bonds                 # martini code for drawing bonds (also highly useful for long atomistic bonds
gmx_dump: @martini/bin/gmxdump      # using GMX 5 with cg_bonds.tcl requires a wrapper around `gmx dump`
martini structure color: false      # martini proteins look best with structure colors (in place of cartoon)
#itp: s02-protein/Protein_A.itp     # structure color requires an itp
backbone color: [None,'black'][1]   # separate color for the BB beads in licorice
cursor color: red                   # protein sidechain color

"""},

'vmd_original_for_reference':{
#####
####
###
##
#
'quick':'view_routine()',
'extensions':[],'tags':[],'imports':['@vmd'],'params':None,
'settings':"""

step: v02-look
video name: video

view mode: ['video','live','snapshot'][0]

#---video settings
max frames: 300
video size: 30
duration: 0.0

#---inputs are autodetected
input gro: system.gro

#---resolution
viewbox: (800,800)
resolution: (2400,1800)
which view: xview
scale: scale by 1.75

#---if Tachyon missing try snapshots
use snapshot: false

#---explicit bonds via cg_bonds.tcl
#cg_bonds: ~/libs/cg_bonds.tcl
#bonder: false

#---special colors for MARTINI proteins
#itp: s02-protein/Protein_A.itp
#backbone color: [None,'black'][1]
#martini structure color: false
#cursor color: red

#---choose a preset "recipe" from vmdmake
recipe collection:|[
	'video aamd atomistic bilayer protein',
	'live aamd atomistic bilayer protein',
	'video cgmd bilayer protein',
	'video cgmd bilayer protein backbone',
	][0]

#---additional highlights
selections:|[
	{'basic_residues':'protein and (resname HIS or resname ARG or resname LYS)',
	'smooth':True,'style':'licorice','goodsell':True},
	][:]

backbone align name: "name CA"

"""},

}