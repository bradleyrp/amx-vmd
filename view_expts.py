{

'vmd_protein':{
#####
####
###
##
#
'quick':"view_routine()",
'extensions':[],'tags':['aamd'],'imports':['@vmd'],'params':None,
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
'extensions':[],'tags':['cgmd','tested_2017.09.20'],'imports':['@vmd'],'params':None,
'settings':"""

USAGE NOTES:|
	this quick script makes a short video of the last gromacs MD part (xtc)
	it is mainly used for reviewing simulations quickly. use omnicalc for videos for publication
	when running locally (i.e. no entry in machine_configuration) make sure to load gromacs to draw bonds

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
cg_bonds: @martini/bin/cg_bonds.tcl # martini code for drawing bonds (also highly useful for long atomistic bonds
gmx_dump: @martini/bin/gmx          # cg_bonds.tcl in @martini/bin takes gmx not gmxdump with the gmx flag
martini structure color: false      # martini proteins look best with structure colors (in place of cartoon)
itp: None                           # structure color requires an itp
backbone color: [None,'black'][1]   # separate color for the BB beads in licorice
cursor color: red                   # protein sidechain color

"""},

'inspect_martini':{
#####
####
###
##
#
'quick':"""view_routine()""",
'extensions':[],'tags':['cgmd','tested_2017.09.20'],'imports':['@vmd'],'params':None,
'settings':"""
USAGE NOTES:|
	this quick script creates a "live" VMD view of a MARTINI simulation
	it saves time by setting up VMD to more easily see a CGMD simulation
step: v99-inspect
view mode: live
input gro: system.gro
viewbox: (800,800)
which view: xview
cg_bonds: @martini/bin/cg_bonds.tcl
gmx_dump: @martini/bin/gmxdump
itp: None # autodetect that last-modified itp file
backbone color: black
martini structure color: true
cursor color: red
recipe collection: live cgmd bilayer protein backbone

"""},

'dextran_pentamer_video':{
#####
####
###
##
#
'quick':"view_routine()",
'extensions':[],'tags':['aamd','tested_2017.09'],'imports':['@vmd'],'params':None,
'settings':"""

USAGE NOTES:|
	not for use with coarse-grained dextran (try vmd_cgmd_bilayer which works fine)
	designed for a brief visualization of the atomistic pentamer dextran
	development note: sometimes the pentamer is off-center

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
backbone align name: "resname AGLC"

#---custom selections (see lib_vmdmake for details)
selections:| [{'basic_residues':'protein and (resname HIS or resname ARG or resname LYS)',
  'smooth':True,'style':'licorice','goodsell':True},
  ][1:0]
#---^^^ customizations turned off if the list is sliced over 1:0"!
#---note: do not use smooth on some items but not others (induces a non-physical asynchronicity)
#---note: SEE lib_vmdmake.py and the vmdmake documentation for more details

DEVNOTES:|
	vmdmake codes need to be updated
	script was modified to the following (quick script method needs updated)
		import os,sys
		sys.path.insert(0,"./runner")
		from makeface import import_remote
		globals().update(**import_remote("inputs/vmd"))
		import amx
		batch.state = amx.state
		batch.settings = amx.settings
		batch.expt = amx.expt
		batch.make_sidestep = amx.make_sidestep
		batch.get_last_frame = amx.get_last_frame
		batch.get_trajectory = amx.get_trajectory
		view_routine()

"""},

}