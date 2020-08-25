#!/usr/bin/env python
#
# Written Aug 14 2020 by J. Pierel 
# Utilities to translate legacy SNANA input files to
# their refactored versions with YAML. Input files are
# for Simulation, LightCurveFit, and SALT2mu/BBC.

import os, sys, yaml, re
import psutil
from   copy		import copy

# definitions
SIM_default_yaml_sections__ = \
							['CONFIG','GENVERSION_LIST','GENOPT_GLOBAL']
SIM_ignore_dict_setup__  = \
						['BATCH_INFO','RANSEED_REPEAT', 'RANSEED_CHANGE', 'RESET_CIDOFF', 
						 'CLEANUP_FLAG', 'GENPREFIX:', 'NGEN_UNIT', 'FORMAT_MASK' ]

SIM_yaml_translation_dict__ = { 'LEGACY':'REFAC' }

SIM_multi_option_list = \
        ['SIMGEN_INFILE_Ia', 'SIMGEN_INFILE_SNIa', 'SIMGEN_INFILE_NONIa' ]

FIT_yaml_translation_dict__ = \
         {'[*]' : '/*/', 
          'LEGACY' : 'REFAC',  
          'APPEND_TABLE_TEXT'   : 'APPEND_TABLE_VARLIST', 
          'FITRES_COMBINE_FILE' : 'APPEND_TABLE_TEXTFILE'
         }

FIT_multi_option_list = [ 'VERSION', 'FITOPT']

BBC_yaml_translation_dict__ = {'XXX*':'blank', 'LEGACY':'REFAC' }

BBC_multi_option_list = []

# ====================================

class _MyDumper(yaml.SafeDumper):
	# HACK: insert blank lines between top-level objects
	# inspired by https://stackoverflow.com/a/44284819/3786245
	def write_line_break(self, data=None):
		super().write_line_break(data)

		if len(self.indents) == 1:
			super().write_line_break()


def _finput_abspath(finput):
	"""
	Helper function that returns absolute path to file
	"""
	finput = finput.strip()
	if not finput.startswith('/') and not finput.startswith('$') and \
	   '/' in finput: finput = '%s/%s'%(cwd,finput)
	return finput

def _has_handle(fpath):
	"""
	Helper function for checking file is safe to open.
	"""

	for proc in psutil.process_iter():
		try:
			for item in proc.open_files():
				if fpath == item.path:
					print(proc.open_files())
					return True
		except Exception:
			pass

	return False

def _open_shared_file(filename,flag="r",max_time=100):
	"""
	Helper function to read a shared file.
	"""

	status = False
	total_time = 0
	while status is False and total_time < max_time:
		if not _has_handle(_finput_abspath(filename)):		  
			f = open(filename, flag)
			status = True
			return f
		else:
			time.sleep(5)
			total_time += 5
	if status is False:
		raise RuntimeError('File %s is opened by another process' %filename)


def _add_keyword_to_dict(current_dict,key,value,legacy_type):
	"""
	Helper function to create nested dictionaries for YAML.
	"""
	kv_list = value.split()
	if len(kv_list)>1 and key not in SIM_ignore_dict_setup__ and\
									 legacy_type=='SIM':
		if key not in current_dict.keys():
			current_dict[key]={}
		if kv_list[1].strip().isnumeric():
			kv_list[1] = int(kv_list[1].strip())
		else:
			try:
				kv_list[1] = float(kv_list[1].strip())
			except:
				kv_list[1] = kv_list[1].strip()
		if kv_list[0].strip() in SIM_multi_option_list:
			current_dict[key][kv_list[0].strip()] = [value.replace(kv_list[0],'').strip()] if\
													len(kv_list) > 2 else\
													[kv_list[1]]
		else:
			current_dict[key][kv_list[0].strip()] = value.replace(kv_list[0],'').strip() if\
													len(kv_list) > 2 else\
													kv_list[1]
	else:
		if value.strip().isnumeric():
			value = int(value.strip())
		else:
			try:
				value = float(value.strip())
			except:
				value = value.strip()
		if key in current_dict.keys():
			if not isinstance(current_dict[key],list):
				current_dict[key] = [current_dict[key]]
		
			current_dict[key].append(value)
		elif any([key in tempList for tempList in [SIM_multi_option_list,
												   FIT_multi_option_list,
												   BBC_multi_option_list]]):
			current_dict[key] = [value]
		else:
			current_dict[key] = value

	return(current_dict)

def _make_yaml_translation(key,value,yaml_key,yaml_value):
	"""
	Helper function that replaces keys between legacy and refactored
	"""

	yaml_key_final = copy(yaml_key)
	yaml_value_final = copy(yaml_value)
	if '*' in yaml_key:
		split_key = yaml_key.split('*')
		split_key_subs = [x for x in split_key if len(x) > 0]
		if not all([substring in key for substring in split_key_subs]):
			actual_key = None
		else:
			if len(split_key[0]) == 0:
				if len(split_key[-1]) == 0:
					actual_key=copy(key)
				else:
					actual_key = key[:key.find(split_key_subs[-1])+len(split_key_subs[-1])]			
			elif len(split_key[-1]) == 0:
				actual_key = key[key.find(split_key_subs[0]):]
			else:	
				actual_key = key[key.find(split_key_subs[0]):key.find(split_key_subs[-1])+len(split_key_subs[-1])]
		if not all([substring in value for substring in split_key_subs]):
			actual_value = None
		else:
			
			if len(split_key[0]) == 0:
				if len(split_key[-1]) == 0:
					actual_value=copy(value)
				else:
					actual_value = value[:value.find(split_key_subs[-1])+len(split_key_subs[-1])]
			elif len(split_key[-1]) == 0:
				actual_value = value[value.find(split_key_subs[0]):]
			else:	
				actual_value = value[value.find(split_key_subs[0]):value.find(split_key_subs[-1])+len(split_key_subs[-1])]

		wildcard_inds_key = [m.start() for m in re.finditer(r'[*]', yaml_key_final)]
		wildcard_inds_value = [m.start() for m in re.finditer(r'[*]', yaml_value_final)]

		yaml_key_final = list(yaml_key_final)
		yaml_value_final = list(yaml_value_final)

		if actual_key == yaml_key.strip('*') or actual_value==yaml_value.strip('*'):
			yaml_key_final = "".join(yaml_key_final)
			yaml_value_final = "".join(yaml_value_final)
			yaml_key_final = yaml_key_final.strip('*')
			yaml_value_final = yaml_value_final.strip('*')
		else:
			for i in range(len(split_key)-1):
				if len(split_key[i]) == 0:
					if actual_key is not None:
						if '*' in yaml_key_final:
							yaml_key_final[wildcard_inds_key[i]]=\
								actual_key[:actual_key.find(split_key[i+1])]
						if '*' in yaml_value_final:
							yaml_value_final[wildcard_inds_value[i]]=\
								actual_key[:actual_key.find(split_key[i+1])]
					elif actual_value is not None:
						if '*' in yaml_key_final:
							yaml_key_final[wildcard_inds_key[i]]=\
								actual_value[:actual_value.find(split_key[i+1])]
						if '*' in yaml_value_final:
							yaml_value_final[wildcard_inds_value[i]]=\
								actual_value[:actual_value.find(split_key[i+1])]
				elif len(split_key[i+1]) == 0:
					if actual_key is not None:
						if '*' in yaml_key_final:
							yaml_key_final[wildcard_inds_key[i]]=\
								actual_key[actual_key.find(split_key[i])+len(split_key[i]):]
						if '*' in yaml_value_final:			
							yaml_value_final[wildcard_inds_key[i]]=\
								actual_key[actual_key.find(split_key[i])+len(split_key[i]):]
					elif actual_value is not None:
						if '*' in yaml_key_final:
							yaml_key_final[wildcard_inds_key[i]]=\
								actual_value[actual_value.find(split_key[i])+len(split_key[i]):]
						if '*' in yaml_value_final:
							yaml_value_final[wildcard_inds_key[i]]=\
								actual_value[actual_value.find(split_key[i])+len(split_key[i]):]
					
				else:
					if actual_key is not None:
						if '*' in yaml_key_final:
							yaml_key_final[wildcard_inds_key[i]]=\
								actual_key[actual_key.find(split_key[i])+len(split_key[i]):actual_key.find(split_key[i+1])]
						if '*' in yaml_value_final:
							yaml_value_final[wildcard_inds_key[i]]=\
							actual_key[actual_key.find(split_key[i])+len(split_key[i]):actual_key.find(split_key[i+1])]
					if actual_value is not None:
						if '*' in yaml_key_final:
							yaml_key_final[wildcard_inds_key[i]]=\
								actual_value[actual_value.find(split_key[i])+len(split_key[i]):actual_value.find(split_key[i+1])]
						if '*' in yaml_value_final:
							yaml_value_final[wildcard_inds_key[i]]=\
								actual_value[actual_value.find(split_key[i])+len(split_key[i]):actual_value.find(split_key[i+1])]

			yaml_key_final = "".join(yaml_key_final)
			yaml_value_final = "".join(yaml_value_final)

	if yaml_key_final in key:
		key = key.replace(yaml_key_final,yaml_value_final)
	if yaml_key_final in value:
		value = value.replace(yaml_key_final,yaml_value_final)
	return(key,value)

def _legacy_snana_sim_input_to_dictionary(basefilename,verbose):
	"""
	Converts a LEGACY sim-input file into a YAML-friendly dictionary.

	Parameters
	----------
	basefilename: str
	The LEGACY sim-input filename

	Returns
	-------
	yaml_dict: dict
	The dictionary to be dumped by YAML
	verbose: bool
	Verbosity flag
	"""

	if not os.path.isfile(basefilename):
		raise ValueError("basefilename cannot be None")
	if verbose:
		print("Load legacy sim input file..",basefilename)
	with _open_shared_file(basefilename) as basefile:
		lines = basefile.readlines()
	
	yaml_dict = {section:{} for section in SIM_default_yaml_sections__}
		
		
	yaml_dict['GENVERSION_LIST']=[]
	inside_genversion_list=False
	for i,line in enumerate(lines):
		if ":" in line and not line.strip().startswith("#"):
			kwline = line.split(":",maxsplit=1)
			kw = kwline[0]
			kv = kwline[1]
			for yaml_key in SIM_yaml_translation_dict__.keys():
				kw,kv=_make_yaml_translation(kw,kv,yaml_key, SIM_yaml_translation_dict__[yaml_key])
			if "#" in kwline[1]:
				kv = kwline[1].split("#")[0].strip()
			if kw == 'GENVERSION':
				yaml_dict['GENVERSION_LIST'].append({'GENVERSION':kv.strip()})
				inside_genversion_list = True
			elif kw == 'ENDLIST_GENVERSION':
				inside_genversion_list = False
			elif inside_genversion_list:
				yaml_dict['GENVERSION_LIST'][-1] = \
						_add_keyword_to_dict(yaml_dict['GENVERSION_LIST'][-1],
											 kw.strip(), kv,'SIM')
					
			elif kw.strip() in SIM_default_yaml_sections__ and kw != 'GENVERSION_LIST':
				yaml_dict = _add_keyword_to_dict(yaml_dict, kw.strip(), kv,'SIM')
			else:
				yaml_dict['CONFIG'] = _add_keyword_to_dict(yaml_dict['CONFIG'],
														   kw.strip(),kv,'SIM')
	
	for top_level in list(yaml_dict.keys()):
		if len(yaml_dict[top_level])==0:
			yaml_dict.pop(top_level)
	return(yaml_dict)

def _legacy_snana_NML_to_dictionary(basefilename,verbose):
	"""
	Converts a LEGACY sim-input file into a YAML-friendly dictionary.

	Parameters
	----------
	basefilename: str
	   The LEGACY lcfit NML filename
	verbose: bool
	  Verbosity Flag

	Returns
	-------
	yaml_dict: dict
	   The dictionary to be dumped by YAML
	"""

	
	if verbose:
		print("Load base fit input file..", basefilename)

	with _open_shared_file(basefilename) as basefile:
		lines = basefile.readlines()

	basekws = []
	nml_lines = []

	snlcinp_fitinp = False
	header = {'CONFIG':{}}
	for i,line in enumerate(lines):
		if '&snlcinp' in line.lower() or '&fitinp' in line.lower(): 
			snlcinp_fitinp = True
		elif '&end' in line.lower():
			snlcinp_fitinp = False
			nml_lines.append(line)
			continue
		if snlcinp_fitinp:
			nml_lines.append(line)
			continue
		

		
		if line.startswith('#'): continue
		if not ':' in line: continue
		key = line.split(':')[0].replace(' ','')
		value = line.split(':')[1].replace('\n','')
		for yaml_key in FIT_yaml_translation_dict__.keys():
			key,value=_make_yaml_translation(key,value,yaml_key,
						FIT_yaml_translation_dict__[yaml_key])

		header['CONFIG'] = _add_keyword_to_dict(header['CONFIG'],key,value,'FIT')

	return(header,nml_lines)


def _legacy_snana_bbc_to_dictionary(basefilename,verbose):
	"""
	Converts a LEGACY sim-input file into a YAML-friendly dictionary.

	Parameters
	----------
	basefilename: str
		The LEGACY sim-input filename

	Returns
	-------
	yaml_dict: dict
		The dictionary to be dumped by YAML
	verbose: bool
	   Verbosity flag
	"""

	if not os.path.isfile(basefilename):
		raise ValueError("basefilename cannot be None")
	if verbose:
		print("Load legacy SALT2mu file..",basefilename)

	with _open_shared_file(basefilename) as basefile:
		lines = basefile.readlines()
	
	yaml_dict = {'CONFIG':{}}
	bbc_lines = []
	inside_genversion_list=False
	for i,line in enumerate(lines):
		if ":" in line and not line.strip().startswith("#"):
			kwline = line.split(":",maxsplit=1)
			kw = kwline[0]
			kv = kwline[1]
			if "#" in kwline[1]:
				kv = kwline[1].split("#")[0].strip()
			for yaml_key in BBC_yaml_translation_dict__.keys():
				kw,kv=_make_yaml_translation(kw,kv,yaml_key,
											 BBC_yaml_translation_dict__[yaml_key])
			yaml_dict['CONFIG'] = _add_keyword_to_dict(yaml_dict['CONFIG'],
													   kw.strip(),kv,'BBC')

		elif '=' in line and not line.strip().startswith('#'):
			bbc_lines.append(line)
	
	for top_level in list(yaml_dict.keys()):
		if len(yaml_dict[top_level])==0:
			yaml_dict.pop(top_level)
	return(yaml_dict,bbc_lines)


# ===================================================
# ===================================================

def SIM_legacy_to_refac(legacy_filename,refactored_filename):
	"""
	Converts a LEGACY SNANA sim-input file into a REFACTORED sim-input file

	Parameters
	----------
	legacy_filename: str
		The filename of the legacy file.
	refactored_filename: str
	   The filename of the refactored file.
	verbose: bool
	   Verbosity flag
	"""

	verbose		= True
	legacy_dict = _legacy_snana_sim_input_to_dictionary(basefilename=legacy_filename,verbose=verbose)
	if verbose:
		print('Outputting refactored YAML file ... %s' % refactored_filename)
		with open(refactored_filename, 'w') as o:
			o.write("# Translated automatically from " \
					"legacy input file %s\n" % legacy_filename )
			yaml.dump(legacy_dict, o, default_flow_style=False,
					  sort_keys=False, Dumper=_MyDumper,width=1000)
			o.write('\n#END_YAML')

#  - - - - - 

def FIT_legacy_to_refac(legacy_filename,refactored_filename):
	"""
	Converts a LEGACY SNANA LCFIT NML file into a REFACTORED LCFIT file

	Parameters
	----------
	legacy_filename: str
		The filename of the legacy file.
	refactored_filename: str
	   The filename of the refactored file.
	verbose: bool
		Verbosity flag
	"""

	verbose = True 
	legacy_header_dict,nml_lines = \
		_legacy_snana_NML_to_dictionary(basefilename=legacy_filename,verbose=verbose)

	if verbose:
		print('Outputting refactored YAML file.. %s'%refactored_filename)

	with open(refactored_filename, 'w') as o :
		o.write('# Automatic translation for legacy LCFIT file %s\n\n' \
				% legacy_filename)
		yaml.dump(legacy_header_dict, o, default_flow_style=False,
				  sort_keys=False, Dumper=_MyDumper,width=1000)

		o.write('\n#END_YAML\n\n')
		for line in nml_lines:
			o.write(line)


def BBC_legacy_to_refac(legacy_filename,refactored_filename ):
	"""
	Converts a LEGACY SNANA BBC file into a REFACTORED BBC file

	Parameters
	----------
	legacy_filename: str
	   The filename of the legacy file.
	refactored_filename: str
	   The filename of the refactored file.
	"""

	verbose = True 
	legacy_dict, bbc_input = \
	_legacy_snana_bbc_to_dictionary(basefilename=legacy_filename,verbose=verbose)
	if verbose:
		print('Outputting refactored YAML file.. %s'%refactored_filename)

	with open(refactored_filename, 'w') as o :
		o.write('# Automatic translation for legacy BBC file %s\n\n' \
				% legacy_filename)
		yaml.dump(legacy_dict, o, default_flow_style=False, 
				  sort_keys=False, Dumper=_MyDumper,width=1000)
		o.write('\n#END_YAML\n\n')
		for line in bbc_input:
			o.write(line)

# END FILE






