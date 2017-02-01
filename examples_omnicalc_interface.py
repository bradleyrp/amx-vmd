#!/usr/bin/python -i

"""
VIDEO-MAKER for molecular simulations
uses vmdmake.VMDWrap on data organized by OMNICALC
designed for atomistic lipid bilayers with proteins
"""

#---plot prep
from plotter import *
from base.store import plotload
import numpy as np
from base.store import plotload,load
sys.path.insert(0,'calcs')
from codes.vmdmake_acme import VMDWrap
from base.gromacs import gmxpaths
import collections,shutil

#---settings
film_folder = 'films'
slice_name = 'survey_short'
group = 'all'
collection = 'normals'

#---prepare directories and get simulation names
out_dn = work.plotdir+'/%s/'%film_folder
if not os.path.isdir(out_dn): os.mkdir(out_dn)
sns = work.load_specs()['collections'][collection]

#---one video per simulation
for sn in sns[:1]:
	#---interface with omnicalc
	work.cursor = (work.c,'xtc')
	gro,traj = [os.path.join(work.postdir,work.slice(sn)[slice_name][group][suf]) 
		for suf in ['gro',work.trajectory_format]]
	tpr = work.get_last_start_structure(sn,part_name='tpr')
	#---reset the cursor in case you want to check work.slice
	work.cursor = (work.c,'xtc')
	if len(work.slice(sn)[slice_name][group]['timeseries'])>300: 
		raise Exception('too many frames. needs development. make a reprocessor!')
	#---find a space to render
	this_dn = os.path.join(out_dn,sn)
	if os.path.isdir(this_dn): raise Exception('refusing to overwrite')
	else: os.mkdir(this_dn)
	outname = 'vid.%s'%re.findall('^.+\/(.+)\.xtc$',traj)[0]
	#---render with VMDWrap
	v = VMDWrap(site=this_dn,gro=gro,xtc=traj,tpr=tpr,res=(1600,1600),last=None,backcolor='white')
	v.do(*'pbc load standard'.split())
	v.trajectory(target=traj)
	v.select(protein='protein and not hydrogen',style='Licorice 0.800000 12.000000 12.000000',
		smooth=True,goodsell=True)
	#---reset the view on just the protein so it is centered (it matters when you do this)
	v.do('ortho','reset','yview')
	v.select(pip2='resname PI2P and not hydrogen',style='licorice',smooth=True,goodsell=True)
	v.select(lipids='not water and not ion and not hydrogen and not protein and not resname PI2P',
		style='licorice',smooth=True,xy=False,goodsell=True,resname_color=True)
	#---try using load then switch to load_dynamic after choosing a good zoom level
	v.command('scale by 0.9')
	v.video()
	v.show(quit=True,clean=False,prompt=False,text=True)
	v.render(name=outname)
