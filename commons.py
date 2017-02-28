
"""
This module contains important constants used by VMDMake.
"""

commons = """
pbc: package require pbc
bonder:|
	source CGBONDSPATH
	cg_bonds -cutoff BOND_CUTOFF_ANGSTROMS -gmx GMXDUMP -tpr TPR
load: mol new GRO
ortho:
	display projection Orthographic
standard:|
	display projection Orthographic
	color Display Background BACKCOLOR
	axes location off
	display resize VIEWX VIEWY
	mol delrep 0 top
dev:|
	mol selection "SELECT"
	mol addrep top
play: animate forward
reset:|
	display resetview
	after 3000
xview:|
	mouse stoprotation
	rotate x to -90
	rotate y by -90
yview:|
	rotate z to 180
	rotate x by -90
zview:|
	mouse stoprotation
	rotate z to 180
backbonename: name CA
align_backbone:|
	set reference [atomselect top "BACKBONENAME" frame 0]
	set compare [atomselect top "BACKBONENAME"]
	set entire [atomselect top "all"]
	set num_steps [molinfo top get numframes]
	for {set frame 0} {$frame < $num_steps} {incr frame} {
		$compare frame $frame
		$entire frame $frame
		set trans_mat [measure fit $compare $reference]
		$entire move $trans_mat
	}
snapshot:|
	set filename SNAPSHOT_FILENAME
	render Tachyon $filename "/usr/local/lib/vmd/tachyon_LINUXAMD64" \\
	-aasamples 12 %s -format TARGA -o %s.tga -res XRES YRES
	exec convert $filename.tga $filename.png
	exec rm $filename.tga $filename
draw_arrow:|
	proc vmd_draw_arrow {mol start end} {
		set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
		graphics $mol cylinder $start $middle radius 0.15
		graphics $mol cone $middle $end radius 0.25
	}
"""

drawing = """
beads: Beads 1.000000 12.000000
vdw: VDW 0.600000 12.000000
licorice: Licorice 0.300000 12.000000 12.000000
licorice_coarse: Licorice 1.300000 12.000000 12.000000
cartoon: NewCartoon 0.300000 10.000000 4.100000 0
cpk_water: CPK 1.000000 0.300000 10.000000 10.000000
lines: Lines 4.000000
points: Points 20.000000
"""

extras = """
update: mol selupdate REP top 1
smooth: mol smoothrep top REP 5
rainbow: mol modcolor REP top ResID
resname_color: mol modcolor REP top ResName
structure_color: mol modcolor REP top Structure
beta_color: mol modcolor REP top Beta
transparent: mol modmaterial REP top Transparent
edgyglass: mol modmaterial REP top EdgyGlass
glass1: mol modmaterial REP top Glass1
goodsell: mol modmaterial REP top Goodsell
hide: mol showrep top REP 0
color_specific: mol modcolor REP top ColorID $color_cursor
glass: mol modmaterial REP top Diffuse
xy:|
	mol showperiodic top REP x
	mol numperiodic top REP 1
	mol showperiodic top REP xy
	mol numperiodic top REP 1
	mol showperiodic top REP xyX
	mol numperiodic top REP 1
	mol showperiodic top REP xyXY
	mol numperiodic top REP 1
pbcz:|
	mol showperiodic top REP z
	mol numperiodic top REP 1
	mol showperiodic top REP zZ
	mol numperiodic top REP 1
"""

defaults = """
res: (1000,1000)
viewbox: (800,800)
backcolor: white
bond_cutoff_angstroms: 20
step: 1
script_fn: script-vmd.tcl
vmd_prepend: "VMDNOCUDA=1"
backbone_select: "name CA"
"""

recipes = """
atomistic bilayer licorice:|
	{'lipids':'not water and not ions and not protein and not hydrogen',
	'smooth':True,'style':'licorice','goodsell':True}
atomistic bilayer lines:|
	{'lipids':'not water and not ions and not protein and not hydrogen',
	'smooth':True,'style':'lines'}
atomistic protein:|
	{'protein_subject':'protein','smooth':True,'style':'cartoon','structure_color':True,'goodsell':True}
atomistic all:|
	{'everything':'all','smooth':True,'style':'cartoon','structure_color':True,'goodsell':True}
waterdots:|
	{'water':'water and name OW','smooth':True,'style':'cpk_water','goodsell':True}
coarse bilayer lines:|
	{'lipids':'not resname W and not resname ION and '+
		'not name BB and not name SC1 and not name SC2 and not name SC3 and not name SC4',
	'smooth':True,'style':'lines'}
coarse bilayer licorice:|
	{'lipids':'not resname W and not resname ION and '+
		'not name BB and not name SC1 and not name SC2 and not name SC3 and not name SC4',
	'smooth':True,'style':'licorice_coarse','goodsell':True}
coarse protein:|
	{'protein_subject':' or '.join(['name %s'%i for i in ['BB','SC1','SC2','SC3','SC4']]),
	'smooth':True,'style':'licorice_coarse','color_specific':True}
coarse protein lines:|
	{'protein_subject':' or '.join(['name %s'%i for i in ['BB','SC1','SC2','SC3','SC4']]),
	'smooth':True,'style':'lines','color_specific':True}
coarse protein backbone:|
	{'protein_subject':' or '.join(['name %s'%i for i in ['BB']]),
	'smooth':True,'style':'licorice_coarse','color_specific':True,'goodsell':True}
coarse protein backbone lines:|
	{'protein_subject':' or '.join(['name %s'%i for i in ['BB']]),
	'smooth':True,'style':'lines'}
"""

recipes_collect = """
video aamd atomistic bilayer protein:['atomistic_bilayer licorice','atomistic_protein']
live aamd atomistic bilayer protein:['atomistic_bilayer lines','atomistic_protein']
video cgmd bilayer protein:['coarse bilayer licorice','coarse protein']
video cgmd bilayer protein backbone:['coarse bilayer licorice','coarse protein backbone']
live cgmd bilayer protein:['coarse bilayer lines','coarse protein']
live cgmd bilayer protein backbone:['coarse bilayer lines','coarse protein backbone lines']
video aamd atomistic all waterdots:['atomistic_bilayer licorice','waterdots']
video aamd generic:['atomistic_bilayer licorice']
"""
