#!/usr/bin/env python

"""
Automatic execution of VMDWrap for the acme-inflected automacs.
Note that alternate runner scripts can be used to interface with other systems like automacs classic, omnicalc
Or better yet you can just write your own VMDWrap-flavored scripts (see examples for examples)
"""

import os,sys,re,json

def create_unbroken_trajectory(vmdstate):
	"""
	Get the last frame from the state and make an unbroken trajectory.
	Note that the typ
	"""
	#---note that typical use of get_last_frame is for chaining/continuation
	#---here we assume that the last state is the one we want, and it still resides in state.json
	get_last_frame(dest=vmdstate.here,source=state.here)
	vmdstate.trajectory_details = get_trajectory(dest=vmdstate.here)
	return vmdstate

def prepare_viewer_scripts(vmdstate):
	"""
	key variables
		show_trajectory : necessary for loading a video
	"""
	kwargs = dict(site=vmdstate.here,gro=vmdstate.trajectory_details['gro'],viewbox=settings.viewbox)
	if settings.get('resolution',None): kwargs['res'] = settings.resolution
	view = VMDWrap(**kwargs)
	cursor_color = settings.get('cursor_color',False)
	if cursor_color: view.set_color_cursor(cursor_color)
	#---custom variables
	custom_variables = settings.get('custom_variables',None)
	if custom_variables:
		assert type(custom_variables)==dict
		view.__dict__.update(**custom_variables)
	#---if dynamics we prepare a trajectory with whole molecules
	#---! currently we have removed show_trajectory from user control
	show_trajectory = settings.view_mode in ['video','live']
	stepskip = 1
	if show_trajectory:
		#---send trajectory files names
		view.__dict__.update(**vmdstate.trajectory_details)
		import warnings
		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			import MDAnalysis
		uni = MDAnalysis.Universe(vmdstate.trajectory_details['gro'],vmdstate.trajectory_details['xtc'])
		max_frames = settings.get('max_frames',None)
		if not max_frames: stepskip = 1
		else: stepskip = max(int(float(len(uni.trajectory))/max_frames+1),1)
		view.__dict__.update(step=stepskip,frames='')
	view.do('load')
	if show_trajectory: view.trajectory(target=vmdstate.trajectory_details['xtc'],step=stepskip)
	view.do('standard')
	backbone_align_name = settings.get('backbone_align_name',False)
	if backbone_align_name:
		view.__dict__['backbonename'] = backbone_align_name
		view.do('align_backbone')
	if settings.get('cg_bonds',False):
		cg_bonds_tcl = os.path.abspath(os.path.expanduser(str(settings.cg_bonds)))
		if not os.path.isfile(cg_bonds_tcl): 
			raise Exception("cannot find cg_bonds.tcl at %s"%cg_bonds_tcl+
				' get this file from MARTINI and set the settings variable called "cg_bonds" to that path')
		view.__dict__['cgbondspath'] = cg_bonds_tcl
		#---! need a protocol for GMXDUMP possibly tempfile or alias
		view.__dict__['gmxdump'] = os.path.abspath(os.path.expanduser(settings.get('gmx_dump','gmxdump')))
		view.do('bonder')
	view.do(settings.which_view)
	if settings.get('scale',False): view.command(settings.scale)
	for select in settings.get('selections',[]): view.select(**select)
	recipe_name = nospaces(settings.recipe_collection)
	if recipe_name not in view.recipes_collect:
		raise Exception('cannot find recipe "%s"'%recipe_name)
	for recipe in view.recipes_collect[nospaces(settings.recipe_collection)]: 
		view.recipe(nospaces(recipe))
	#---apply special MARTINI protein colors
	if settings.get('martini_structure_color',False): 
		martini_style = 'backbone'
		#---! previously allowed another option in some cases: martini_style = 'all'
		view.martini_secondary_structure(itp=settings.itp,
			style=martini_style,back_color=settings.backbone_color)
		#---!method for polymers needs to be fixed up a bit
		if False:
			view.martini_secondary_structure(itp=settings.itp,style=martini_style,
				back_color=settings.backbone_color,mers=settings.get('mers',1))
	#---custom scripts directly from the wordspace
	if settings.get('custom_script',False):
		script = settings.custom_script
		assert type(script)==list,'custom script must be a list'
		for line in script: eval(line)
	#---handle snapshots
	snap_dn_abs = 'snapshots'
	#---stage 3: make snapshots if this is a video
	if settings.view_mode == 'video' and not settings.get('video_snaps',False): 
		settings.video_snaps = view.video(dn=snap_dn_abs,snapshot=settings.get('use_snapshot',False))
	elif settings.view_mode == 'video': 
		view.video_dn,view.snapshot = settings.video_snaps,settings.get('use_snapshot',False)
	text_mode = settings.view_mode in ['video','snapshot']
	#---finished
	vmdstate.view = view
	return vmdstate

def run_vmd(vmdstate):
	"""
	"""
	view = vmdstate.view
	text_mode = settings.view_mode in ['video','snapshot']
	view.show(quit=text_mode,text=text_mode)
	vmdstate.already_ran = True
	return vmdstate

def render_video(vmdstate):
	"""
	"""
	view = vmdstate.view
	if settings.view_mode == 'video': 
		view.render(name=settings.video_name,size=settings.video_size,
			duration=settings.get('duration',0),webm=settings.webm)
	return vmdstate

def view_routine():
	"""
	Manage incoming video requests.
	This function performs the four steps required to make images with VMD.
	It is designed to rerender videos without rerendering snapshots if desired.
	"""

	for key in ['state','settings','expt']:
		if key not in globals(): raise Exception('%s must be exported to this module'%key)

	#---prepare the vmdstate
	vmdstate = DotDict(**{'status':'started'})
	state_fn = os.path.join(settings.step,'vmdstate.json')
	#---load the state if it exists
	if os.path.isfile(state_fn): 
		with open(state_fn) as fp: vmdstate = DotDict(**json.load(fp))
	if 'here' not in vmdstate: make_sidestep(settings.step)
	vmdstate.here = os.path.join(settings.step,'')

	try:
		if not vmdstate.get('already_ran',False): vmdstate = create_unbroken_trajectory(vmdstate)
		vmdstate = prepare_viewer_scripts(vmdstate)
		if not vmdstate.get('already_ran',False): vmdstate = run_vmd(vmdstate)
		vmdstate = render_video(vmdstate)
	#---on failure or interrupt we save the state
	except (Exception,KeyboardInterrupt) as e: 
		print('[VMD] exception: %s'%str(e))
		vmdstate['status'] = 'error'
		#---! tracebackstolen from acme.stopper
		import sys,traceback
		exc_type,exc_obj,exc_tb = sys.exc_info()
		#---write the traceback
		tag = '[TRACEBACK] '
		tracetext = tag+re.sub(r'\n','\n%s'%tag,str(''.join(traceback.format_tb(exc_tb)).strip()))
		print(tracetext)
	print('[VMD] writing state to %s'%state_fn)
	with open(state_fn,'w') as fp: 
		if 'view' in vmdstate: vmdstate.pop('view')
		fp.write(json.dumps(vmdstate))
