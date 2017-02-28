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
cg_bonds: @cg_bonds                 # martini code for drawing bonds (also highly useful for long atomistic bonds
gmx_dump: @martini/bin/gmxdump      # using GMX 5 with cg_bonds.tcl requires a wrapper around `gmx dump`
martini structure color: false      # martini proteins look best with structure colors (in place of cartoon)
#itp: s02-protein/Protein_A.itp     # structure color requires an itp
backbone color: [None,'black'][1]   # separate color for the BB beads in licorice
cursor color: red                   # protein sidechain color

"""},

'inspect_martini':{
#####
####
###
##
#
'quick':'view_routine()',
'extensions':[],'tags':[],'imports':['@vmd'],'params':None,
'settings':"""
USAGE NOTE:|
	currently needs to work on snapshots
	also we have to wait for state.json to be written, so maybe write it before production starts?
step: v99-inspect
view mode: live
input gro: system.gro
viewbox: (800,800)
which view: xview
cg_bonds: ~/libs/cg_bonds.tcl
gmx_dump: @martini/bin/gmxdump
itp: s02-adhere/Protein.itp
backbone color: black
martini structure color: true
cursor color: red
recipe collection: live cgmd bilayer protein backbone

"""},

}