#!/usr/bin/env python

import os,sys,re,tempfile,json

__all__ = ['martiniview']

def martiniview(*args,**kwargs):
	"""
	Quickly view a martini structure.
	!Needs more fine-tuning.
	!Currently deprecated?
	"""
	#---! hardcoded-path
	cg_bonds_fn = '~/libs/cg_bonds.tcl'
	gmxdump_fn = 'inputs/martini/bin/gmxdump'
	#---! first we check for state.json
	if os.path.isfile('state.json'):
		sys.path.insert(0,'')
		import amx
		last = amx.get_last_gmx_call('mdrun')
		gro,tpr,xtc = [os.path.join(amx.state.here,last['flags'][i]) for i in ['-c','-s','-x']]
	else:
		def parse(ext):
			matches = [i for i in args if os.path.splitext(i)[1]==('.'+ext)]
			if len(matches)>1: raise Exception('too many matches to extension %s in the command line'%ext)
			if not matches: return None
			return matches[0]
		gro = kwargs.get('gro',parse('gro'))
		tpr = kwargs.get('tpr',parse('tpr'))
		xtc = kwargs.get('xtc',parse('xtc'))
	if not all([gro,tpr]): raise Exception('missing one of gro (%s) tpr (%s)'%(gro,tpr))
	cmds = ['mol new %s'%gro]
	if xtc: cmds += ['mol addfile %s'%xtc]
	cmds += ['source %s'%cg_bonds_fn]
	cmds += ['cg_bonds -gmx %s -tpr %s'%(gmxdump_fn,tpr)]
	tf = tempfile.NamedTemporaryFile(delete=False)
	tf.write(('\n'.join(cmds)+'\n').encode())
	tf.close()
	cmd = 'vmd -e %s'%tf.name
	print('[STATUS] running "%s"'%cmd)
	os.system(cmd)

def lammps_view_basic():
	"""
	Mimics view_routine for lammps.
	""" 

	import amx
	#---manually sending functions from automacs to vmdmake
	#---previously checked for these in globals
	for key in ['state','settings','expt',
		'make_sidestep','gmx_get_last_frame','gmx_get_trajectory']:
		globals()[key] = getattr(amx,key)

	#---prepare the vmdstate
	vmdstate = DotDict(**{'status':'started'})
	state_fn = os.path.join(settings.step,'vmdstate.json')
	#---load the state if it exists
	if os.path.isfile(state_fn): 
		with open(state_fn) as fp: vmdstate = DotDict(**json.load(fp))
	if 'here' not in vmdstate: make_sidestep(settings.step)
	vmdstate.here = os.path.join(settings.step,'')
	try:
		if not vmdstate.get('already_ran',False):
			#---if the user supplies a custom trajectory then we use that
			if settings.get('trajectory_details',False): 
				vmdstate.trajectory_details = settings['trajectory_details']
				#---paths checked here and should be relative to root.
				for key,val in vmdstate.trajectory_details.items():
					if not os.path.isfile(val): 
						raise Exception('incoming trajectory_details file missing: %s'%val)
					else: vmdstate.trajectory_details[key] = os.path.abspath(val)
			else:
				#---automatically fix PBCs on the last part if the user has not supplied a custom file
				vmdstate = create_unbroken_trajectory(vmdstate)
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
