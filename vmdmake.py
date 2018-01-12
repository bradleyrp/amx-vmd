#!/usr/bin/python

import os,sys,json,time,re,subprocess,tempfile,glob

#---(acme) magic for local import when you run from elsewhere
sys.path.insert(0,os.path.dirname(os.path.relpath(os.path.abspath(__file__),os.getcwd())))

from commons import *

#---! development note: add cog/com from http://www.ks.uiuc.edu/Research/vmd/vmd-1.7.1/ug/node181.html

#---yamlb requires no spaces in keys
def nospaces(text): return re.sub(' ','_',text)

class VMDWrap:
	#---hard-coded vmd colors
	vmd_colors = dict([(i,ii) for ii,i in enumerate(
		('blue red gray orange yellow tan silver green white pink cyan purple lime '+
		'mauve ochre iceblue black yellow2 yellow3 green2 green3 cyan2 cyan3 blue2 blue3 violet '+
		'magenta magenta2 red2 red3 orange2 orange3').split())])

	def yamlparse(self,text):
		"""Try yamlparse (comes from yamlb in acme) first otherwise treat it like python."""
		#---we ignore JSON errors so we can use python code and also use the "tabbed" style in yamlb
		try: out = yamlparse(text,style='tabbed',ignore_json=True)
		except: out = eval(text)
		return out

	def __init__(self,name='vmd',**kwargs):
		"""
		Interface with VMD to make videos and execute common rendering functions.
		Need to ensure viewbox corresponds to resolution.
		"""
		global commons,defaults,extras,drawing
		self.commons = self.yamlparse(commons)
		self.defaults = self.yamlparse(defaults)
		self.extras = self.yamlparse(extras)
		self.drawing = self.yamlparse(drawing)
		self.recipes = self.yamlparse(recipes)
		self.recipes_collect = self.yamlparse(recipes_collect)

		#---add any tcl scripts in this folder to the commons
		for fn in glob.glob(os.path.join(os.path.dirname(__file__),'*.tcl')):
			tcl_lib = re.search('^(.+)\.tcl$',os.path.basename(fn)).group(1)
			assert tcl_lib not in self.commons,'already found "%s" in the commons'%tcl_lib
			self.commons[tcl_lib] = open(fn).read()

		#---load defaults
		self.__dict__.update(**self.defaults)
		#---all kwargs are saved to the dictionary as substitutions
		self.__dict__.update(**kwargs)
		self.start_time = time.time()
		#---incoming resolution overrides the default viewbox specified in commons.py
		if 'res' in self.__dict__.keys(): 
			self.res = tuple(self.res)
			self.viewbox = self.res
		self.viewx,self.viewy = self.viewbox

		#---handle paths where the 'site' is the top level folder
		self.cwd = kwargs.get('site',os.getcwd())
		if not os.path.isdir(self.cwd): raise Exception('requested site does not exist: %s'%self.cwd)

		#---counters
		self.script = []
		self.selections = {}
		self.selection_order = []
		self.set_color_cursor('black')

	def __setitem__(self,key,value):
		"""
		Use the class like a dictionary for updating substitutions in the commons.
		"""
		self.__dict__[key] = value

	def do(self,*names):
		"""
		Use a set of prescribed code blocks defined in commons.
		"""
		for name in names:
			if name not in self.commons: 
				raise Exception('cannot find item "%s" in the commons/constants'%name)
			chunk = self.commons[name].strip('\n')
			for key,val in self.__dict__.items(): 
				chunk = re.sub('(\s|\")(%s)(\s|$|\")'%key.upper(),r"\g<1>%s\g<3>"%str(val),chunk)
			self.script.append(chunk)

	def trajectory(self,target,first=0,last=-1,step=1):
		"""
		Load a trajectory.
		Most functions are handled in the commons.py, however it makes more sense to do this manually.
		"""
		self.command('animate delete beg 0 end 0 skip 0 0')
		if not os.path.isfile(target): raise Exception('cannot find file %s'%target)
		self.command("mol addfile %s first %d last %d step %d waitfor all"%(target,first,last,step))

	def select(self,**kwargs):
		"""
		Create a new selection with string substitutions.
		"""
		style = kwargs.pop('style') if 'style' in kwargs else None
		extras = []
		for key in list(kwargs.keys()):
			if type(kwargs[key]) == bool:
				val = kwargs.pop(key)
				if val: extras.append(key)
		dynamics = kwargs.pop('dynamics') if 'dynamics' in kwargs else None
		for key in kwargs:
			selstring = str(kwargs[key])
			for sel in self.selections:
				selstring = re.sub(sel.upper(),'(%s)'%self.selections[sel],selstring)
			self.script += ["mol selection \"(%s)\""%selstring,"mol addrep top"]		
			self.script += ['set %s [atomselect top "%s"]'%(key,selstring)]
			self.selection_order.append(key)
			self.selections[key] = kwargs[key]
		if style != None:
			for key in kwargs:
				style = self.drawing[style] if style in self.drawing else style
				self.script += ['mol modstyle %d top %s'%(self.selection_order.index(key),style)]
		for extra in extras:
			if extra in self.extras and extra:
				for key in kwargs:
					instruct = str(self.extras[extra])
					instruct = re.sub('REP','%d'%self.selection_order.index(key),instruct)
					self.script += [instruct]
			else: raise Exception('missing extra setting: %s'%extra)			

	def set_color_cursor(self,color):
		"""
		The $color_cursor variable in tcl is set to black (16) by the standard package but can be changed
		to use coloring by a particular ColorID, on the fly.
		"""
		color = {'grey':'gray'}.get(color,color)
		if color in self.vmd_colors: self.script.append('set color_cursor %d'%self.vmd_colors[color])
		else: self.script.append('set color_cursor %s'%color)

	def command(self,text): 
		"""
		Add a raw command to the script if necessary.
		"""

		self.script.append(text)

	def write(self):
		"""
		Write the script.
		"""

		with open(self.cwd+'/'+self.script_fn,'w') as fp: 
			for line in self.script: fp.write(line+'\n')

	def recipe(self,recipe,**kwargs):
		"""
		Add from a set of predefined selections.
		"""

		recipe = dict(self.recipes[recipe])
		recipe.update(**kwargs)
		self.select(**recipe)

	def show(self,**kwargs):
		"""
		Run VMD.
		The pipe option will feed commands directly to VMD for systems with a wonky setup.
		!option to remove everything when done (clean)?
		!dispdev text option (text)?
		"""

		quit = kwargs.get('quit',True)
		pipe = kwargs.get('pipe',False)
		#---text mode if quit on complete
		text = kwargs.get('text',quit)

		#---quit after the script is done
		if quit: self.script += ['quit']
		self.write()
		#---! arrest execution if you on light running dark: import ipdb;ipdb.set_trace()
		#---open VMD for the user		
		if quit == False:
			#---if we are not quitting then we run vmd directly
			#---note this might fail if the weird segfault that motivated "not run_script" happens
			os.system('cd %s && %s vmd %s -e %s'%(self.cwd,self.vmd_prepend,
				'-dispdev text' if text else '',self.script_fn))
		#---quit after running
		else: 
			#---feed the commands directly to VMD to preclude the weird segfault errors
			if not pipe:
				text = '\n'.join(self.script)
				proc = subprocess.Popen('vmd -dispdev text',stdin=subprocess.PIPE,
					shell=True,stdout=sys.stdout,stderr=sys.stdout,
					cwd=self.cwd)
				proc.communicate(input=str(text).encode())
			else:
				#---call the VMD script (note that this fails sometimes on systems with a poor VMD setup)
				subprocess.check_call('%s vmd %s -e %s'%(self.vmd_prepend,'-dispdev text' 
					if text else '',self.script_fn),shell=True,cwd=self.cwd)
		return True

	def video(self,traces=None,snapshot=False,pause=False,cwd='',dn=None):

		"""
		Prepare to render a video by writing an add-on script.
		"""

		#---ensure that the viewbox and the resolution have the same proportions
		if float(self.viewbox[0])/self.viewbox[1]!=float(self.res[0])/self.res[1]:
			self.viewbox = tuple([int(round(float(max(self.viewbox))/max(self.res)*i)) for i in self.res])
			self.script.append('display resize %d %d'%self.viewbox)
		#---add a "trace" or tail to some elements
		if traces != None:
			trace_commands = []
			if type(traces) != list: raise Exception('traces must be a list or None')
			for key in traces:
				trace_commands.append('mol drawframes top %d $lag:$lead'%self.selection_order.index(key))
		else: trace_commands = []
		#---generic header for making a video
		video_lines = [
			'set num [molinfo top get numframes]',
			'puts "STATUS number of frames is $num"',
			'for {set i 0} {$i < $num} {incr i} {',
			'animate goto $i',
			'puts "STATUS advanced to frame $i/$num"',
			'set lag [expr $i-15]',
			'set lead [expr $i]']
		video_lines.extend(trace_commands)
		#---all videos have snapshots rendered in a temporary folder
		if not dn: self.video_dn = tempfile.mkdtemp(dir=self.cwd)
		else:
			self.video_dn = dn
			dropspot = os.path.join(self.cwd,self.video_dn)
			if not os.path.exists(dropspot): os.mkdir(dropspot)
			######### assert glob.glob(os.path.join(dropspot,'*')) == []
		dropspot = os.path.join(os.path.basename(self.video_dn),'')
		#---snapshot is the cheap/fast method otherwise we render
		self.snapshot = snapshot
		if not self.snapshot:
			#---get the generic path for tachyon
			video_lines.extend([
				"set archexe \"\"",
				"switch [vmdinfo arch] {",
				"	WIN64 -","	WIN32 { set archexe \".exe\" }",
				"}","package require exectool 1.2","set tachyonexe [format \"tachyon%s\" $archexe];",
				"set tachyoncmd [::ExecTool::find -interactive -description \"Tachyon Ray Tracer\" "+
				"-path [file join $env(VMDDIR) \"tachyon_[vmdinfo arch]$archexe\"] $tachyonexe]"])
			video_lines.extend([
				'set filename '+dropspot+'snap.[format "%04d" $i]',
				'render TachyonInternal $filename.tga' if False else
				'render Tachyon $filename $tachyoncmd '+\
				'-aasamples 12 %s -format TARGA -o %s.tga -res '+'%d %d'%self.res,
				'exec convert $filename.tga $filename.png',
				'exec rm $filename.tga $filename',
				'}'])
		else:
			video_lines.extend([
				'set filename '+dropspot+'snap.[format "%04d" $i]',
				'render snapshot $filename.ppm',
				'}'])
		#---write script-video.tcl which will be sourced by script-vmd.tcl
		with open(self.cwd+'/script-video.tcl','w') as fp:
			for line in video_lines: fp.write(line+'\n')
		if not pause: self.script.append('source %sscript-video.tcl'%os.path.join(cwd,''))
		return self.video_dn

	def render(self,name='video',rate=1.0,codec=None,size=50.0,duration=0.0,webm=False):
		"""
		Render a video from snapshots.
		"""

		def is_installed(command):
			std,err = subprocess.Popen(command.split(),shell=True,executable='/bin/bash',
				stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
			return not re.search('command not found','\n'.join([str(std),str(err)]))

		encoders = ['ffmpeg','avconv']
		use_encoder = next((i for i in encoders if is_installed(i)),None)
		if not use_encoder: raise Exception('cannot find any encoders in %r'%encoders)

		#---when using ffmpeg confirm that the codec is available 
		if not codec and use_encoder=='ffmpeg':
			#---priority codecs
			ffmpeg_codec_priorities = [('mpeg4','mp4'),('mpeg2video','mpeg')]
			proc = subprocess.Popen('ffmpeg -codecs',shell=True,
				stdout=subprocess.PIPE,stdin=subprocess.PIPE)
			std,err = proc.communicate()
			if sys.version_info>=(3,0): std = std.decode()


			"""
			note that ffmpeg -codecs reports some header material, a line with some dashes, and then a 
			several-column list of codecs where the second and third columns are the code and the name
			"""
			try: start_lineno = next(ii for ii,i in enumerate(std.splitlines()) if re.match('^\s*\-+\s*$',i))
			except: raise Exception('cannot locate a ---- in the `ffmpeg -codecs` output')
			regex_col1 = '^\s*[^\s]+\s+([^\s]+)'
			codecs = [re.match(regex_col1,j).group(1) for j in 
				std.splitlines()[start_lineno:] if re.match(regex_col1,j)]
			codecs_found = [(n,ext) for n,ext in ffmpeg_codec_priorities if n in codecs]
			if not codecs_found: raise Exception('codec list is %s\n'%str(codecs)+
				'but we cannot find any requested codecs %s'%(str(ffmpeg_codec_priorities)))
			else: codec = codecs_found[0][0]
			video_suffix = dict(ffmpeg_codec_priorities)[codec]

		assert not (size and use_encoder!='ffmpeg')
		two_pass = '2pass' if size else '1pass'

		#---index commands by encoder,self.snapshot
		dropspot = os.path.join(os.path.basename(self.video_dn),'')
		#---bitrate for 2pass is currently somewhat conservative and undershoots the desired size
		nframes = float(len(glob.glob(os.path.join(self.cwd,self.video_dn,'*.png'))))
		frate = 24.0
		if not size: bitrate = None 
		elif duration==0.0: 
			try: bitrate = float(size)*8192/(float(nframes)/frate)
			except: raise Exception('invalid bitrate: size=%s, nframes=%s, frate=%s'%(size,nframes,frate))
		#---duration overrides rates
		else: 
			bitrate = float(size)*8192/(float(duration))
			rate = float(duration)/(nframes/frate)
		commands = {
			('ffmpeg','snapshot','1pass'):r"ffmpeg -i "+dropspot+"snap.%04d.ppm "+"-vcodec %s -q 0"%codec+
				"-filter:v 'setpts=%.1f*PTS' "%rate+name+"."+video_suffix,
			('ffmpeg','tachyon','1pass'):r"ffmpeg -i "+dropspot+"/snap.%04d.png "+"-vcodec %s -q 0 "%codec+
				"-filter:v 'setpts=%.1f*PTS' "%rate+name+"."+video_suffix,
			#---codec options for avconv are currently hardcoded
			('avconv','snapshot','1pass'):r"avconv -r %.1f -i "%(rate*4)+dropspot+
				"/snap.%04d.png "+self.cwd+'/'+name+".mp4",
			('avconv','tachyon','1pass'):r"avconv -r %.1f -i "%(rate*4)+dropspot+
				"/snap.%04d.png "+self.cwd+'/'+name+".mp4",
			#---! currently the 2pass method is hard-coded
			('ffmpeg','tachyon','2pass'):r"ffmpeg -y -i "+dropspot+"snap.%04d.png "+
				"-framerate %d "%frate+"-filter:v 'setpts=%.1f*PTS' "%round(rate,3)+
				"-c:v libx264 -b:v %dk -pass 1 -f mp4 /dev/null && "%int(bitrate)+
				"ffmpeg -y -i "+dropspot+"snap.%04d.png -framerate "+"%d "%frate+
				"-c:v libx264 -preset medium -b:v %dk -pass 2 "%int(bitrate)+
				"-filter:v 'setpts=%.1f*PTS' "%round(rate,3)+name+".mp4",
			}
		command_key = (use_encoder,'snapshot' if self.snapshot else 'tachyon',two_pass)
		assert rate==1.0 or command_key==('ffmpeg','tachyon','2pass'),'rates must be 1.0 for 2pass'
		if command_key not in commands: 
			raise Exception('cannot parse the following command options for the encoder: "%r"'%command_key)
		subprocess.check_call(commands[command_key],cwd=self.cwd,shell=True)
		print('[NOTE] ffmpeg command was: "%s"'%commands[command_key])
		if command_key==('ffmpeg','tachyon','2pass'):
			print('[NOTE] cleaning up ffmpeg2pass files')
			for fn in glob.glob(self.cwd+'/ffmpeg2pass*'): os.remove(fn)
		#---option for a webm conversion when complete
		if webm:
			subprocess.check_call("ffmpeg -i %s -c:v libvpx -crf 10 -b:v 1M -c:a libvorbis %s"%
				(name+"."+video_suffix,name+".webm"),shell=True,cwd=self.cwd)
			
	def martini_secondary_structure(self,itp,style='all',name='protein_subject',back_color=None,mers=1):
		"""
		Color a MARTINI protein by secondary structure.
		"""
		import re
		#---set the jet colorscale from colorscales.tcl
		#---! needs added to vmdmake self.do('colorscales')
		#---! needs added to vmdmake self.script.append('colorscale_jet')
		#---backbone forces secondary structure coloring
		if style=='backbone': back_color = None
		#---select secondary structure colors
		coloring = {'C':'gray','H':'red','3':'pink','1':'pink','S':'blue','2':'pink','T':'yellow'}
		#---! HACKED
		coloring = {'C':0.25,'H':1.0,'3':1.0,'1':1.0,'S':0.0,'2':1.0,'T':0.55}
		#coloring = {'C':0.0,'H':1.0,'3':0.0,'1':0.0,'S':0.0,'2':0.0,'T':0.0}
		#---hard-coded bead numbers generated from lib_helices.CoarseGrainedModel
		bead_counts = {'ILE':2,'GLN':2,'GLY':1,'PHE':4,'GLU':2,'CYS':2,'HIS':4,'SER':2,
			'LYS':3,'PRO':2,'HISH':4,'ASN':2,'VAL':2,'THR':2,'ASP':2,'ASP0':2,'LYS0':3,
			'ARG0':3,'GLU0':2,'ALA':1,'MET':2,'LEU':2,'ARG':3,'TYR':4,'TRP':5}
		aa_codes3 = ['TRP','TYR','PHE','HIS','ARG','LYS','CYS','ASP','GLU',
			'ILE','LEU','MET','ASN','PRO','HYP','GLN','SER','THR','VAL','ALA','GLY']
		aa_codes1 = 'WYFHRKCDEILMNPOQSTVAG'
		#---autodetect the ITP by running a python find command and sorting by mtime
		if not itp:
			#---simple way to get the last file
			try: itp = sorted([os.path.join(rn,fn) for rn,dn,fns in os.walk('.') 
				for fn in fns if re.match('^.+\.itp$',fn)],key=lambda x:os.path.getmtime(x))[-1]
			except: raise Exception('failed to autodetect the itp file. try setting that flag in the experiment')
		if not os.path.isfile(itp):
			raise Exception('cannot locate ITP %s'%itp)
		with open(itp) as fp: text = fp.read()
		regex_ss = '^;\s*(S|s)econdary structure\n;\s*([A-Z0-9]+)\s'
		ss = re.search('^;\s*(?:S|s)econdary (?:S|s)tructure.*?\n;\s*(.*?)$',text,re.M+re.DOTALL).group(1)
		seq = re.search('^;\s*(?:S|s)equence.*?\n;\s*(.*?)$',text,re.M+re.DOTALL).group(1)
		mapping = dict([(i,float(ii)) for ii,i in enumerate(['C','H','3','1','S','2','T'])])
		obs_bead_counts = [bead_counts[dict(zip(aa_codes1,aa_codes3))[i]] for i in seq]
		#betas = [(self.vmd_colors[back_color] if i==0 and back_color else self.vmd_colors[coloring[s]]) 
		#	for ii,s in enumerate(ss) for i in range({'all':obs_bead_counts[ii],'backbone':1}[style])]
		betas = [coloring[s] for ii,s in enumerate(ss) 
			for i in range({'all':obs_bead_counts[ii],'backbone':1}[style])]
		#---if the ITP is a single monomer and there is more than one, we use the "mers" flag
		if mers>1: betas = betas*mers
		betas_string = "{%s}"%' '.join(['%.2f'%i for i in betas])
		self.script.append('set beta_%s %s'%(name,betas_string))
		self.script.append('$%s set beta $beta_%s'%(name,name))
		self.script.append('mol modcolor %d top Beta'%self.selection_order.index('protein_subject'))
