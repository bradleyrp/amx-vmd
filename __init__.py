#!/usr/bin/env python

import os,sys

from . import vmdmake
from . import batch
from . vmdmake import VMDWrap
from . batch import view_routine
from . batch_original import view_routine_original

#---search for the acme runner and if found send the state to vmdmake
#---note that the search for config.py mimics the __init__.py that pairs with importer.py for top-level
#---...acme-style modules, however we treat this module (vmd) like a standard module instead
#---...as a result, we manually import some important features of a standard automacs/acme run
#---...the use of a standard import avoids using "from amx import *" to get e.g. the state

acme_root = False
root_dns = [os.path.abspath(os.path.join(__file__,*('..' for j in range(i+1)))) for i in range(5)]
try: acme_root = next(f for f in root_dns if os.path.isfile(os.path.join(f,'config.py')))
except: print('[WARNING] running vmdmake without acme')
config = {}

#---if we find an ACME config file we import several useful functions
#---note that we may find an omnicalc config.py when using vmdmake inside omnicalc
#---...which case is apparent when the acme key is missing in which case we skip the imports for automacs
if acme_root and config.get('acme',False):

	with open(os.path.join(acme_root,'config.py')) as fp: config = eval(fp.read())
	#---connect to runner
	sys.path.insert(0,os.path.abspath(os.path.dirname(config['acme'])))

	#---amx dependencies (new acme-version functions, so we import them iff acme)
	import amx

	vmdmake.get_last_frame = batch.get_last_frame = amx.get_last_frame
	vmdmake.get_trajectory = batch.get_trajectory = amx.get_trajectory
	vmdmake.make_sidestep = batch.make_sidestep = amx.make_sidestep

	#---acme dependencies
	from runner import state,settings,expt
	from runner.datapack import yamlb,DotDict

	vmdmake.state = batch.state = state
	vmdmake.settings = batch.settings = settings
	vmdmake.expt = batch.expt = expt

#---codes for operation without ACME
else: 
	try: from acme_supplement import DotDict,yamlb
	except: print('[WARNING] no DotDict or yamlb')

#---use a copy of the ACME simplified YaML parser
vmdmake.DotDict = batch.DotDict = DotDict
vmdmake.yamlparse = yamlb

#---we have to export VMDWrap (instead of importing it) since it depends on yamlb
batch.VMDWrap = VMDWrap
batch.nospaces = vmdmake.nospaces
