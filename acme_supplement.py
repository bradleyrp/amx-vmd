#!/usr/bin/env python

import os,sys,shutil,json,glob,tempfile,subprocess,re,ast,pprint,time
from makeface_supplement import abspath,strip_builtins,import_remote,str_or_list,fab

class DotDict(dict):
	"""Use dots to access dictionary items."""
	def __init__(self,*args,**kwargs):
		#---special sauce that maps the dictionary self to its own attributes list
		self.__dict__ = self
		#---maintain a list of protected keywords used by the dict object
		self.__dict__['_protect'] = list(dict.__dict__.keys())
		#---load the dictionary with incoming args and kwargs
		for key,val in args: self.__dict__[key] = val
		self.update(kwargs)
	def __setattr__(self,key,val):
		if key in self.get('_protect',[]): 
			raise Exception('cannot use dict/protected attribute %s'%key)
		#---special sauce that allows you to set new attributes
		super(dict,self).__setattr__(key,val)
	def __repr__(self): 
		#---hide the underscore attributes
		return str(dict([(i,j) for i,j in self.items() if not i.startswith('_')]))
	def __getattr__(self,key):
		"""
		The above code is sufficient for a standard DotDict.
		The following allows fallback lookups.
		Development note: add protection that asserts that the query function q is callable
		Note that attribute lookups here throw an error if no match (because there is no default),
		hence you should use the query function directly for "get" (with default) functionality.
		"""
		if key in self: return self[key]
		#---failed attributes are looked up with the query function if available
		elif 'q' in self: return self.q(key)
		else: raise KeyError('cannot find key %s'%str(key))

def jsonify(text): 
	"""
	Convert python to JSON by replacing single quotes with double quotes and stripping trailing commas
	Note that this checker might be oversensitive if additional JSON errors with the nested dictionary
	Used before SafeDictHook to check for redundant keys. We use a placeholder in block text because JSON cannot 
	process it and the point of this function is to use SafeDictHook to prevent redundant keys.
	"""
	#---remove comments because they might screw up the JSON
	text = re.sub(r'([\"]{3}.*?[\"]{3})','"REMOVED_BLOCK_COMMENT"',text,flags=re.M+re.DOTALL)
	#---note that this fails if you use hashes inside of dictionary values
	text = re.sub(r'(#.*?)\n','',text,flags=re.M+re.DOTALL)
	#---strip trailing commas because they violate JSON rules
	text = re.sub(r",[ \t\r\n]*}","}",text.replace("'","\""))
	#---fix the case on all booleans
	text = re.sub("True","true",text)
	text = re.sub("false","false",text)
	#---! rpb is worried that this is a hack
	return text

class SafeDictHook(dict):
	"""Hook for json.loads object_pairs_hook to catch repeated keys."""
	def __init__(self,*args,**kwargs):
		self.__class__ == dict
		if len(args)>1: raise Exception('development failure')
		keys = [i[0] for i in args[0]]
		if len(keys)!=len(set(keys)): raise Exception('repeated keys (or JSON error) in %s'%str(keys))
		self.update(*args,**kwargs)

def check_repeated_keys(text,verbose=False):
	"""Confirm that dict literals pass through non-redundant json checker."""
	extra_msg = "either fix the repeated keys or check for JSON problems."
	text_json = jsonify(text)
	try: _ = json.loads(text_json,object_pairs_hook=SafeDictHook)
	except Exception as e: 
		print('[ERROR] found repeated keys (or JSON encoding error). '+extra_msg)
		if verbose: 
			print('[ERROR] for your reference, the jsonify string is: \"\"\"\n%s\n\"\"\"'%text_json)
			print('[ERROR] exception is %s'%e)
		return False
	return True

def checker(fn):
	"""
	Check a file for evaluation from the command line.
	! Less generic name.
	"""
	if not os.path.isfile(fn): raise Exception('cannot find file %s'%fn)
	with open(fn) as fp: text = fp.read()
	print('[NOTE] parsing with python')
	result = eval(text)
	print('[NOTE] parsing with jsonify')
	result = jsonify(text)
	print('[NOTE] parsing with check_repeated_keys')
	check_repeated_keys(text,verbose=True)

def yamlb(text,style=None,ignore_json=False):
	"""
	Basic YAML parser.
	Development note: missing colons are hard to troubleshoot. Predict them?
	Development note: doesn't prevent errors with multiple keys in a dictionary!
	"""
	unpacked,compacts = {},{}
	str_types = [str,unicode] if sys.version_info<(3,0) else [str]
	#---evaluate code blocks first
	regex_block_standard = r"^\s*([^\n]*?)\s*(?:\s*:\s*\|)\s*([^\n]*?)\n(\s+)(.*?)\n(?!\3)"
	regex_block_tabbed = r"^\s*([^\n]*?)\s*(?:\s*:\s*\|)\s*\n(.*?)\n(?!\t)"
	if style == 'tabbed': regex_block = regex_block_tabbed
	else: regex_block = regex_block_standard
	regex_line = r"^\s*(.*?)\s*(?:\s*:\s*)\s*(.+)$"
	#---strip comments first 
	text = re.sub("\s*#.*?$",'',text,flags=re.M)
	while True:
		blockoff = re.search(regex_block,text,re.M+re.DOTALL)
		if not blockoff: break
		if style == 'tabbed': key,block = blockoff.groups()[:2]
		else: 
			#---collect the key, indentation for replacement, and value
			key,indent,block = blockoff.group(1),blockoff.group(3),''.join(blockoff.groups()[1:])
		#---alternate style does multiline blocks with a single tab character
		#---! who uses this? only vmdmake? might be worth dropping
		if style == 'tabbed': compact = re.sub("(\n\t)",r'\n',block.lstrip('\t'),re.M)
		#---remove indentations and newlines (default)
		else: compact = re.sub('\n','',re.sub(indent,'',block))
		unpacked[re.sub(' ','_',key)] = compact
		#---remove the block
		text = re.sub(re.escape(text[slice(*blockoff.span())]),'',text)
	while True:
		line = re.search(regex_line,text,re.M)
		if not line: break
		key,val = line.groups()
		if key in unpacked and not ignore_json:
			raise Exception('\n[ERROR] key is repeated in the settings: "%s"'%key)
		unpacked[re.sub(' ','_',key)] = val
		text = re.sub(re.escape(text[slice(*line.span())]),'',text)
	#---evaluate rules to process the results
	for key,val_raw in unpacked.items():
		#---store according to evaluation rules
		try: val = eval(val_raw)
		except SyntaxError as e:
			#---use of the explicit dict constructor catches repeated keywords
			if e.msg=='keyword argument repeated': raise Exception('repeated keys in: %s'%val_raw)
			else: val = val_raw
		except: val = val_raw
		if type(val)==list: result = val
		elif type(val)==dict:
			if not ignore_json and not check_repeated_keys(val_raw):
				raise Exception('repeated keys in: "%s"'%val_raw)
			result = val
		elif type(val) in str_types:
			if re.match('^(T|t)rue$',val): result = True
			elif re.match('^(F|f)alse$',val): result = False
			#---! may be redundant with the eval command above
			elif re.match('^[0-9]+$',val): result = int(val)
			elif re.match('^[0-9]*\.[0-9]*$',val): result = float(val)
			else: result = val
		else: result = val
		unpacked[key] = result
	return unpacked
