#!/usr/bin/env python

import re,os,tempfile

__all__ = ['martiniview']

def martiniview(*args,**kwargs):
	"""
	Quickly view a martini structure.
	!Needs more fine-tuning.
	"""
	#---! hardcoded-path
	cg_bonds_fn = '~/libs/cg_bonds.tcl'
	gmxdump_fn = 'inputs/martini/bin/gmxdump'
	def parse(ext):
		matches = [i for i in args if os.path.splitext(i)[1]==('.'+ext)]
		if len(matches)>1: raise Exception('too many matches to extension %s in the command line'%ext)
		if not matches: return None
		return matches[0]
	gro = kwargs.get('gro',parse('gro'))
	tpr = kwargs.get('tpr',parse('tpr'))
	xtc = kwargs.get('xtc',parse('xtc'))
	if not all([xtc,gro,tpr]): raise Exception('missing one of gro (%s) tpr (%s)'%(gro,tpr))
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