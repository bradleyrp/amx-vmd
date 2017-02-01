#!/usr/bin/env python

"""
Combined automacs (classic, original) and omnicalc execution for making standard videos.
See `runner.py` for the updated, acme-inflected version.
Or see the examples for some alternate ways of using VMDWrap.
"""

def view_routine_original(**kwargs):
	"""
	Function which interprets incoming video requests for VMDWrap and then executes them.
	Note that this function prepares a local viewer script for each video so users can tune/rerender.
	"""

	raise Exception('this needs back-ported for automacs (classic) and omnicalc')

	#---distinguish incoming paths based on amx/omni calls
	caller = kwargs.pop('caller','amx')
	assert caller in ['omni','amx'],'caller must be omni or amx'
	#---if automacs calls, the wordspace is fully populated
	wordspace = kwargs.get('wordspace',None)
	#---if omnicalc calls, we populate the wordspace
	if not wordspace:
		if caller=='amx': raise Exception('must pass wordspace if running from amx')
		wordspace = WordSpace()
		wordspace['trajectory_details'] = kwargs.pop('trajectory_details',{})
	#---fold settings into the wordspace
	wordspace.update(**yamlparse(kwargs['settings']))
	snap_dn_abs = kwargs.pop('snap_dn_abs','snapshots')
	video_name = kwargs.pop('video_name',wordspace.get('video_name','video'))
	#---amx uses a local wordspace for iterative rerendering
	if caller=='amx': wordspace_out = wordspace.step+'/'+'wordspace.json'
	elif caller=='omni': 
		film_root = os.path.abspath(os.path.join(snap_dn_abs,'..'))
		wordspace_out = os.path.join(film_root,'%s.wordspace.json'%video_name)
		#---if the wordspace file exists we might be doing a rerender
		if os.path.isfile(wordspace_out):
			print('[NOTE] found wordspace for this video: "%s"'%wordspace_out)
			incoming_wordspace = json.load(open(wordspace_out))
			wordspace.update(**incoming_wordspace)
			#---current settings override the previous wordspace
			#---previous wordspace really just tracks the development state for rerendering
			wordspace.update(**yamlparse(kwargs['settings']))

	#---recapitulate script-viewer here with amx/omni switching
	try:
		if wordspace.view_mode=='video': wordspace.show_trajectory = True
		#---stage 1: create trajectory with unbroken molecules
		if not wordspace.get('under_development',False):
			if caller=='amx':
				wordspace['last_step'],wordspace['last_part'] = detect_last()
				#---! previously used "start" function from automacs
				try: start(wordspace.step)
				except: os.mkdir(wordspace.step)
				#---always get the last frame for viewing
				get_last_frame(tpr=True,cpt=False,ndx=True,top=False,itp=False)
				wordspace.trajectory_details = remove_box_jumps(step=wordspace.step,
					last_step=wordspace.last_step,last_part=wordspace.last_part)
			else:
				if 'trajectory_details' not in wordspace:
					raise Exception('view_routine called with caller="amx" '+
						'but we need trajectory_details')
			wordspace.under_development = 1
			write_wordspace(wordspace,outfile=wordspace_out)
		#---stage 2: prepare the viewer scripts
		if caller=='amx': site = wordspace.step
		elif caller=='omni': site = film_root
		view = VMDWrap(site=site,
			gro='system-input.gro' if caller=='amx' else wordspace.trajectory_details['gro'],
			res=wordspace.resolution,viewbox=wordspace.viewbox)
		if wordspace['cursor_color']: view.set_color_cursor(wordspace.cursor_color)
		#---custom variables
		if wordspace['custom_variables']:
			assert type(wordspace.custom_variables)==dict
			view.__dict__.update(**wordspace.custom_variables)
		#---if dynamics we prepare a trajectory with whole molecules
		if wordspace.show_trajectory:
			#---send trajectory files names
			view.__dict__.update(**wordspace.trajectory_details)
			#---use MDAnalysis to get the right frame counts
			import MDAnalysis
			uni = MDAnalysis.Universe(
				(wordspace.step+'/' if caller=='amx' else '')+view.gro,
				(wordspace.step+'/' if caller=='amx' else '')+view.xtc)
			stepskip = max(int(float(len(uni.trajectory))/wordspace.max_frames+1),1)
			view.__dict__.update(step=stepskip,frames='')
		view.do('load_dynamic' if wordspace.show_trajectory else 'load','standard')
		if wordspace['backbone_align_name']:
			view.__dict__['backbonename'] = wordspace.backbone_align_name
			view.do('align_backbone')
		if wordspace['cg_bonds'] or wordspace.bonder:
			cg_bonds_tcl = os.path.abspath(os.path.expanduser(str(wordspace['cg_bonds'])))
			if not os.path.isfile(cg_bonds_tcl): 
				raise Exception('\n'+''.join(['[ERROR] %s\n'%i for i in 
					['cannot find cg_bonds.tcl','get this file from MARTINI',
					'and set the settings variable called "cg_bonds" to that path']]))
			view.__dict__['cgbondspath'] = cg_bonds_tcl
			#---! need a protocol for GMXDUMP possibly tempfile or alias
			view.__dict__['gmxdump'] = 'gmxdump'
		if wordspace.bonder: view.do('bonder')
		view.do(wordspace.which_view)
		if wordspace['scale']: view.command(wordspace.scale)
		for select in wordspace.selections: view.select(**select)
		for recipe in view.recipes_collect[nospaces(wordspace.recipe_collection)]: 
			view.recipe(nospaces(recipe))
		#---apply special MARTINI protein colors
		if wordspace['martini_structure_color']: 
			if wordspace['recipe_collection']=='video cgmd bilayer protein backbone':
				martini_style = 'backbone'
			else: martini_style = 'all'
			view.martini_secondary_structure(itp=wordspace.itp,
				style=martini_style,back_color=wordspace.backbone_color)
			view.martini_secondary_structure(itp=wordspace.itp,style=martini_style,
				back_color=wordspace.backbone_color,mers=wordspace.get('mers',1))
		#---custom scripts directly from the wordspace
		if wordspace['custom_script']:
			script = wordspace.custom_script
			assert type(script)==list,'custom script must be a list'
			for line in script: eval(line)
		#---stage 3: make snapshots if this is a video
		if wordspace.view_mode == 'video' and not wordspace['video_snaps']: 
			wordspace.video_snaps = view.video(dn=snap_dn_abs,snapshot=wordspace['use_snapshot'])
			wordspace.under_development = 2
			write_wordspace(wordspace,outfile=wordspace_out)
		elif wordspace.view_mode == 'video': 
			view.video_dn,view.snapshot = wordspace.video_snaps,wordspace['use_snapshot']
		#---stage 3: run show and create snapshots if making a video
		if wordspace.under_development < 3:
			text_mode = wordspace.view_mode in ['video','snapshot']
			view.show(quit=text_mode,text=text_mode)
			wordspace.under_development = 3
			write_wordspace(wordspace,outfile=wordspace_out)
		#---stage 4: render the video
		if wordspace.view_mode == 'video': 
			view.render(name=video_name,
				size=wordspace.video_size,duration=wordspace.get('duration',0))
		wordspace.under_development = 4 if wordspace.view_mode == 'video' else 2
	except KeyboardInterrupt as e: exception_handler(e,wordspace,all=True)
	except Exception as e: exception_handler(e,wordspace,all=True)

def remove_box_jumps_DEPRECATED():
	"""
	Convert the trajectory to reassemble broken molecules.
	Requires items frm state.last_mdrun_output.
	Note that this is customized for vmdmake but it could be generalized and added to automacs.py.
	"""
	last_tpr = state.last_mdrun_output['-s']
	last_xtc = state.last_mdrun_output['-x']
	last_partno = int(re.match('^md\.part([0-9]{4})',os.path.basename(last_xtc)).group(1))
	trjconv_call = state.gmxpaths['trjconv']
	cmd = trjconv_call+' -f %s -s %s -o %s -pbc mol'%(
		os.path.join('../',last_xtc),last_tpr,
		'md.part%04d.pbcmol.xtc'%last_partno)
	bash(cmd,log='trjconv-pbcmol',cwd=state.here,inpipe='0\n')
	cmd = trjconv_call+' -f %s -s %s -o %s -pbc mol'%(
		'system-input.gro',last_tpr,'system-input.pbcmol.gro')
	bash(cmd,log='trjconv-pbcmol-gro',cwd=state.here,inpipe='0\n')
	gro = 'system-input.gro'
	tpr = 'system-input.tpr'
	xtc = 'md.part%04d.pbcmol.xtc'%last_partno
	return {'xtc':xtc,'gro':gro,'tpr':tpr}
