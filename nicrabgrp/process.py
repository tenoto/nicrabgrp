__author__  = 'Teruaki Enoto'
__version__ = '0.03'
__date__    = '2018 November 26'
"""
HISTORY
2018-11-27 transfered from giantradiopulse to nicrabgrp library
2018-10-24 created by T.Enoto 
"""

import os 
import sys 
import shutil
import glob 
import yaml
import numpy as np 
import pandas as pd 
import astropy.io.fits as fits 

import nicrabgrp.radio
import nicrabgrp.xray
import nicrabgrp.const

class ObservationUnit():
	def __init__(self,row,param,outdir):
		print("-- ObservationUnit {} is generated.".format(row['dataid']))
		self.param = {}
		for keyword in param:
			self.param[keyword] = param[keyword]
		for keyword in row.index:
			self.param[keyword] = row[keyword]	
		self.outdir = outdir 

	def make_suboutdir(self):
		print('--- %s method: %s ---' % (self.param['dataid'],sys._getframe().f_code.co_name))
		self.suboutdir = '%s/%s' % (self.outdir,self.param['dataid'])
		self.param['suboutdir'] = self.suboutdir
		self.param['suboutdir_radio'] = '%s/radio' % self.suboutdir		
		self.param['suboutdir_xray'] = '%s/xray' % self.suboutdir				
		if not os.path.exists(self.param['suboutdir']):
			os.makedirs(self.suboutdir)
		if not os.path.exists(self.param['suboutdir_radio']):
			os.makedirs(self.param['suboutdir_radio'])
		if not os.path.exists(self.param['suboutdir_xray']):
			os.makedirs(self.param['suboutdir_xray'])

	def show_parameters(self):
		print('--- %s method: %s ---' % (self.param['dataid'],sys._getframe().f_code.co_name))
		print(self.param)

	def set_gti_file(self):
		print('--- %s method: %s ---' % (self.param['dataid'],sys._getframe().f_code.co_name))
		search_filename = '%s/%s_*_DE???NEW?_GTI.txt' % (self.param['PATH_TO_RADIO_DATADIR'],self.param['dataid'])
		if len(glob.glob(search_filename)) == 0:
			sys.stderr.write('file {} not found.\n'.format(search_filename))
			return -1
		elif len(glob.glob(search_filename)) > 1:
			sys.stderr.write('file {} has multiple candidates.\n'.format(search_filename))
			return -1
		else:
			self.radio_gtifile_txt = glob.glob(search_filename)[0]
			sys.stdout.write('gti file {} is set.\n'.format(self.radio_gtifile_txt))
			return 0 

	def set_grplist_file(self,grptype='MP'):
		print('--- %s method: %s ---' % (self.param['dataid'],sys._getframe().f_code.co_name))
		search_filename = '%s/%s_*_%sGRPlistDE???SNge*PHpm?_updated.txt' % (self.param['PATH_TO_RADIO_DATADIR'],self.param['dataid'],grptype)
		#search_filename = '%s/%s_*_%sGRPlistDE???SNge*_updated.txt' % (self.param['PATH_TO_RADIO_DATADIR'],self.param['dataid'],grptype)
		if len(glob.glob(search_filename)) == 0:
			sys.stderr.write('file {} not found.\n'.format(search_filename))
			return -1
		elif len(glob.glob(search_filename)) > 1:
			sys.stderr.write('file {} has multiple candidates.\n'.format(search_filename))
			return -1
		else:
			if grptype == 'MP':
				self.mpgrplist_txt = glob.glob(search_filename)[0]
				sys.stdout.write('mpgrplist file {} is set.\n'.format(self.mpgrplist_txt))
			elif grptype == 'IP':
				self.ipgrplist_txt = glob.glob(search_filename)[0]
				sys.stdout.write('ipgrplist file {} is set.\n'.format(self.ipgrplist_txt))				
			else:
				sys.stderr.write('invalid grptype [should be MP or IP]\n')				
				return -1 
			return 0 	

	def convert_radiogti_txt2fits(self):
		print('--- %s method: %s ---' % (self.param['dataid'],sys._getframe().f_code.co_name))
		self.radio_gti = nicrabgrp.radio.GootTimeIntervalTextFile(self.radio_gtifile_txt)
		if 'radio_gti_fitsfile' in self.param:
			sys.stdout.write('radio_gti_fitsfile already exists. skip.\n')
		else:
			self.param['radio_gti_fitsfile'] = self.radio_gti.writeAsFitsFormat()
			shutil.move(self.param['radio_gti_fitsfile'],self.param['suboutdir_radio'])
			self.param['radio_gti_fitsfile'] = '%s/%s' % (self.param['suboutdir_radio'],self.param['radio_gti_fitsfile'])

	def convert_radiompgrp_txt2fits(self):
		print('--- %s method: %s ---' % (self.param['dataid'],sys._getframe().f_code.co_name))
		self.radio_mpgrp = nicrabgrp.radio.GiantRadioPulseTextFile(self.mpgrplist_txt)
		if 'radio_mpgrp_fitsfile' in self.param:
			sys.stdout.write('radio_mpgrp_fitsfile already exists. skip.\n')
		else:					
			self.param['radio_mpgrp_fitsfile'] = self.radio_mpgrp.writeAsFitsFormat()
			shutil.move(self.param['radio_mpgrp_fitsfile'],self.param['suboutdir_radio'])
			self.param['radio_mpgrp_fitsfile'] = '%s/%s' % (self.param['suboutdir_radio'],self.param['radio_mpgrp_fitsfile'])

	def convert_radioipgrp_txt2fits(self):
		print('--- %s method: %s ---' % (self.param['dataid'],sys._getframe().f_code.co_name))
		self.radio_ipgrp = nicrabgrp.radio.GiantRadioPulseTextFile(self.ipgrplist_txt)
		if 'radio_ipgrp_fitsfile' in self.param:
			sys.stdout.write('radio_ipgrp_fitsfile already exists. skip.\n')
		else:		
			self.param['radio_ipgrp_fitsfile'] = self.radio_ipgrp.writeAsFitsFormat()
			shutil.move(self.param['radio_ipgrp_fitsfile'],self.param['suboutdir_radio'])		
			self.param['radio_ipgrp_fitsfile'] = '%s/%s' % (self.param['suboutdir_radio'],self.param['radio_ipgrp_fitsfile'])

	def set_met_epoch(self):
		cmd  = 'rm -rf tmp_nitimeconv.out;'
		cmd += 'nitimeconv.py %d -f mjd -s tt > tmp_nitimeconv.out' % self.param['MJD']
		print(cmd);os.system(cmd)
		with open('tmp_nitimeconv.out') as f:
			for line in f:
				cols = line.split()
				if len(cols) == 0:
					continue 
				if cols[0] == 'NICER':
					met_epoch_day = float(cols[7])
		self.param['met_epoch'] = met_epoch_day + float(self.param['tJPLms'])*1e-3
		sys.stdout.write('MET epoch: {:.16f}\n'.format(self.param['met_epoch']))
		cmd  = 'rm -rf tmp_nitimeconv.out;'
		print(cmd);os.system(cmd)

	def set_nicer_filenames(self):
		path_to_nicer_obsid_candidates = glob.glob('%s/*/%s' % (self.param['PATH_TO_XRAY_DATADIR'],self.param['nicerobsid']))
		if len(path_to_nicer_obsid_candidates) == 0:
			sys.stdout.write('no nicer obsid\n')
			quit()
		self.param['path_to_nicer_obsid'] = path_to_nicer_obsid_candidates[0]
		self.param['nievt_cl_original'] = '%s/xti/event_cl/ni%s_0mpu7_cl.evt.gz' % (self.param['path_to_nicer_obsid'],self.param['nicerobsid'])
		self.param['nievt_ufa_original'] = '%s/xti/event_cl/ni%s_0mpu7_ufa.evt.gz' % (self.param['path_to_nicer_obsid'],self.param['nicerobsid'])		
		self.param['niorbfile'] = '%s/auxil/ni%s.orb.gz' % (self.param['path_to_nicer_obsid'],self.param['nicerobsid'])

		self.param['nievt_cl_bary']  = '%s/ni%s_0mpu7_cl_nibary.evt' % (self.param['suboutdir_xray'],self.param['nicerobsid'])
		self.param['nilog_cl_bary']  = '%s/ni%s_0mpu7_cl_nibary.log' % (self.param['suboutdir_xray'],self.param['nicerobsid'])		
		self.param['nievt_cl_bary_phase'] = '%s/ni%s_0mpu7_cl_nibary_phase.evt'% (self.param['suboutdir_xray'], self.param['nicerobsid'])
		self.param['nievt_cl_bary_phase_grpflag'] = '%s/ni%s_0mpu7_cl_nibary_phase_grpflag.evt'% (self.param['suboutdir_xray'], self.param['nicerobsid'])		
		self.param['nievt_cl_bary_phase_grpflag_esel'] = '%s/ni%s_0mpu7_cl_nibary_phase_grpflag_%sto%skeV.evt'% (self.param['suboutdir_xray'], self.param['nicerobsid'],
			str(self.param['NICER_ENERGY_MIN_KEV']).replace('.','p'),str(self.param['NICER_ENERGY_MAX_KEV']).replace('.','p'))
		self.param['niscript_process'] = '%s/ni%s_process.sh' % (self.param['suboutdir_xray'],self.param['nicerobsid'])
		self.param['niscript_addflag'] = '%s/ni%s_addflag.sh' % (self.param['suboutdir_xray'],self.param['nicerobsid'])

	def make_nicer_process_script(self):
		dump  = '#!/bin/sh -f\n\n'

		dump += '# --------- Barycentric correction --------- \n'
		dump += 'barycorr \\\n'
		dump += 'infile=%s \\\n'  % self.param['nievt_cl_original']
		dump += 'outfile=%s \\\n' % self.param['nievt_cl_bary']		
		dump += 'ra=%.6f dec=%.6f \\\n' % (self.param['CRAB_PULSAR_RA_J2000'],self.param['CRAB_PULSAR_DEC_J2000'])
		dump += 'orbitfiles=%s \\\n' % self.param['niorbfile']
		dump += 'refframe=ICRS ephem=JPLEPH.%s \\\n' % self.param['PLANETARY_EPHEMERIS'].replace('DE','')
		dump += '>& %s\n\n\n' % self.param['nilog_cl_bary']

		dump += '# ---------  add BARY_TIME column --------- \n'
		dump += 'fcalc clobber=yes infile=%s[\'EVENTS\'] ' % self.param['nievt_cl_bary']
		#dump += 'fcalc clobber=yes infile=%s ' % self.param['nievt_cl_bary']
		dump += 'outfile=%s clname="BARY_TIME" ' % self.param['nievt_cl_bary']
		dump += 'expr="#MJDREFI + #MJDREFF + TIME/86400.0" \n\n'

		dump += 'fcalc clobber=yes infile=%s[\'GTI\'] ' % self.param['nievt_cl_bary']
		#dump += 'fcalc clobber=yes infile=%s ' % self.param['nievt_cl_bary']
		dump += 'outfile=%s clname="BARY_START" ' % self.param['nievt_cl_bary']
		dump += 'expr="#MJDREFI + #MJDREFF + START/86400.0" \n\n'

		dump += 'fcalc clobber=yes infile=%s[\'GTI\'] ' % self.param['nievt_cl_bary']
		#dump += 'fcalc clobber=yes infile=%s ' % self.param['nievt_cl_bary']		
		dump += 'outfile=%s clname="BARY_STOP" ' % self.param['nievt_cl_bary']
		dump += 'expr="#MJDREFI + #MJDREFF + STOP/86400.0" \n\n'

		dump += '# ---------  add PULSE_PHASE column --------- \n'
		# This calculation of F2 comes from the C code at the end of the
		# explanatory notes for the Jodrell ephemeris
		# f2 = 2.0*p1*p1/(p0*p0*p0)	
		p0 = float(self.param['P0ms'])*1e-3
		p1 = float(self.param['Pdot0'])
		self.param['f2'] = 2.0*p1*p1/(p0*p0*p0)			
		#dump += 'rm -f %s\\\n' %  self.param['nievt_cl_bary_phase'] 
		dump += '%s \\\n' % shutil.which("faddphase_nu.py")
		dump += '%s \\\n' % self.param['nievt_cl_bary']  # infits             Input event fits file.
		dump += '%.7f \\\n' % self.param['met_epoch']    # epoch              Folding epoch in a unit of TIME.
		dump += '%.13f \\\n' % float(self.param['nu0'])  # nu                 Folding frequency nu (Hz)
		dump += '--nudot=%6.5fe-15 \\\n' % float(self.param['nudot0_15'])
		dump += '--nu2dot=%.6e \\\n' % self.param['f2']		
		dump += '--outfits=%s \\\n' % self.param['nievt_cl_bary_phase'] 

		with open(self.param['niscript_process'],'w') as f:
			f.write(dump)	
		cmd = 'chmod +x %s' % self.param['niscript_process']
		print(cmd);os.system(cmd)

	def run_nicer_process_script(self):
		print('--- %s method: %s ---' % (self.param['dataid'],sys._getframe().f_code.co_name))		
		print('...runnning command...\n')
		cmd = './%s' % self.param['niscript_process']
		print(cmd);os.system(cmd)

	def add_grpflag_to_xrayevents(self):
		print('--- %s method: %s ---' % (self.param['dataid'],sys._getframe().f_code.co_name))		

		out_column_name_mpgrp = 'MPGRP_SN%.1f' % (self.param['SNthr'])
		out_column_name_ipgrp = 'IPGRP_SN%.1f' % (self.param['SNthr'])
		xrayevt = nicrabgrp.xray.XrayEventFile(self.param['nievt_cl_bary_phase'])
		xrayevt.add_grpflag(
			self.param['radio_mpgrp_fitsfile'],
			self.param['radio_ipgrp_fitsfile'],
			self.param['SNthr'],
			self.param['nievt_cl_bary_phase_grpflag'],
			out_column_name_mpgrp,out_column_name_ipgrp
			)

		pi_min = int(self.param['NICER_ENERGY_MIN_KEV'] * nicrabgrp.const.NICER_KEV2PI)
		pi_max = int(self.param['NICER_ENERGY_MAX_KEV'] * nicrabgrp.const.NICER_KEV2PI)
		cmd = 'fselect %s %s "(PI >= %d) && (PI <= %d)"' % (self.param['nievt_cl_bary_phase_grpflag'],
			self.param['nievt_cl_bary_phase_grpflag_esel'],pi_min,pi_max)
		print(cmd);os.system(cmd)

	def plot_pulse_profile(self):
		print('--- %s method: %s ---' % (self.param['dataid'],sys._getframe().f_code.co_name))

		self.param['nipls_cl_bary_phase_grpflag_esel'] = self.param['nievt_cl_bary_phase_grpflag_esel'].replace('.evt','_pls.fits')
		xrayevtlist = nicrabgrp.xray.XrayEventList([self.param['nievt_cl_bary_phase_grpflag_esel']])
		xrayevtlist.generate_profile_fitsfile(self.param['nipls_cl_bary_phase_grpflag_esel'])

		for nphase in nicrabgrp.xray.NPHASE_LIST_DEFAULT:
			outpdf = '%s_n%d.pdf' % (self.param['nipls_cl_bary_phase_grpflag_esel'].replace('.fits',''),nphase)
			title  = '%s, Phase N%d, ' % (self.param['dataid'],nphase)
			title += '{:.1e} X-rays, {:.1e} MP-GRPs ({:.3f}%)'.format(xrayevtlist.numofevt_all,xrayevtlist.numofevt_mpgrp,xrayevtlist.fraction_ipgrp*100.0)
			xrayprofile = nicrabgrp.xray.XrayProfile(self.param['nipls_cl_bary_phase_grpflag_esel'])
			xrayprofile.plot_profile_fitsfile(nphase=nphase,outpdf=outpdf,xmin=0.0,xmax=2.0,ymin=None,ymax=None,title=title)
			outpdf_zoom = '%s_n%d_zoom.pdf' % (self.param['nipls_cl_bary_phase_grpflag_esel'].replace('.fits',''),nphase)			
			xrayprofile.plot_profile_fitsfile(nphase=nphase,outpdf=outpdf_zoom,xmin=0.94,xmax=1.06,ymin=None,ymax=None,title=title)

		enhancement, significance = xrayprofile.get_enhancement_significance(target_bins_list=self.param['GRP_SEARCH_TARGET_BINS_LIST'],nphase=self.param['GRP_SEARCH_NPHASE'])
		self.param['enhancement_N%d'] = enhancement
		self.param['significance_N%d'] = enhancement		
		print(enhancement, significance)

	def write_parameter_yamlfile(self):
		print('--- %s method: %s ---' % (self.param['dataid'],sys._getframe().f_code.co_name))
		outyamlfile = '%s/%s_setup.yaml' % (self.param['suboutdir'],self.param['dataid'])
		f = open(outyamlfile,'w')
		f.write(yaml.dump(self.param, default_flow_style=False))
		f.close()

	def reload_parameter_yamlfile(self,inputyamlfile):
		print('--- %s method: %s ---' % (self.param['dataid'],sys._getframe().f_code.co_name))
		self.inputyamlfile = inputyamlfile

		if not os.path.exists(self.inputyamlfile):
			raise FileNotFoundError("{} not found.".format(self.inputyamlfile))
		try:
			self.param = yaml.load(open(self.inputyamlfile))
		except OSError as e:
			raise 
		print("setup yaml file {} is successfully loaded.".format(self.inputyamlfile))

class ProcessManager():
	""" 
	:param file_path: path to a file to setup yaml file.
	"""

	def __init__(self,file_path,outdir='out/crabgrp'):
		print('--- ProcessManager ---')	
		# set output directory.
		self.outdir = outdir 
		#if not os.path.exists(self.outdir):
		#	os.makedirs(self.outdir)

		# read yaml file.
		self.file_path = file_path		
		if not os.path.exists(self.file_path):
			raise FileNotFoundError("setup yaml file {} not found.".format(self.file_path))
		try:
			self.df = pd.read_csv(self.file_path,header=1,delim_whitespace=True)
		except OSError as e:
			raise 
		print("setup yaml file {} is successfully loaded.".format(self.file_path))
		self.param = yaml.load(open(self.file_path))

	def read_ephemeris_file(self):
		print('--- ProcessManager method: %s ---' % sys._getframe().f_code.co_name)
		self.df = pd.read_csv(self.param['CRAB_PULSAR_EPHEMERIS_FILE'],
			delim_whitespace=True,header=0)

		self.observationunit_list = []
		for index, row in self.df.iterrows():
			if row['proc_flag'] != True:
				sys.stdout.write('-- ObservationUnit {} is skipped.\n'.format(row['dataid']))
				continue 
			self.observationunit_list.append(ObservationUnit(row,self.param,self.outdir))			
			search_filename = '%s/%s/%s_setup.yaml' % (self.outdir,row['dataid'],row['dataid'])
			if len(glob.glob(search_filename)) > 0:
				self.observationunit_list[-1].reload_parameter_yamlfile(glob.glob(search_filename)[0])			
			self.observationunit_list[-1].show_parameters()		

	def convert_radiofiles(self):
		print('--- ProcessManager method: %s ---' % sys._getframe().f_code.co_name)		

		# set individual observations 
		self.read_ephemeris_file()	

		for obs in self.observationunit_list:
			obs.make_suboutdir()
			obs.set_gti_file()		
			obs.set_grplist_file(grptype='MP')
			obs.set_grplist_file(grptype='IP')
			obs.convert_radiogti_txt2fits()
			obs.convert_radiompgrp_txt2fits()
			obs.convert_radioipgrp_txt2fits()
			obs.write_parameter_yamlfile()

	def prepare_xrayfiles(self):	
		print('--- ProcessManager method: %s ---' % sys._getframe().f_code.co_name)		

		# set individual observations 
		self.read_ephemeris_file()	

		for obs in self.observationunit_list:
			obs.make_suboutdir()		
			obs.set_nicer_filenames()
			obs.set_met_epoch()
			obs.make_nicer_process_script()
			obs.run_nicer_process_script()
			obs.write_parameter_yamlfile()

	def add_grpflag_to_xrayfiles(self):
		print('--- ProcessManager method: %s ---' % sys._getframe().f_code.co_name)		

		# set individual observations 
		self.read_ephemeris_file()	

		for obs in self.observationunit_list:
			obs.add_grpflag_to_xrayevents()
			obs.write_parameter_yamlfile()

	def plot_individual_profiles(self):
		print('--- ProcessManager method: %s ---' % sys._getframe().f_code.co_name)		

		# set individual observations 
		self.read_ephemeris_file()	

		for obs in self.observationunit_list:
			obs.show_parameters()
			obs.plot_pulse_profile()
			obs.write_parameter_yamlfile()

	def accumulated_significance(self,indir,outdir):
		print('--- ProcessManager method: %s ---' % sys._getframe().f_code.co_name)		

		self.df = pd.read_csv(self.param['CRAB_PULSAR_EPHEMERIS_FILE'],
			delim_whitespace=True,header=0)

		self.observationunit_list = []
		for index, row in self.df.iterrows():
			if row['add_flag'] != True:
				sys.stdout.write('-- ObservationUnit {} is skipped.\n'.format(row['dataid']))
				continue 
			self.observationunit_list.append(ObservationUnit(row,self.param,indir))			
			search_filename = '%s/%s/%s_setup.yaml' % (indir,row['dataid'],row['dataid'])
			if len(glob.glob(search_filename)) > 0:
				self.observationunit_list[-1].reload_parameter_yamlfile(glob.glob(search_filename)[0])			
			#self.observationunit_list[-1].show_parameters()		

		if outdir != '' and not os.path.exists(outdir):
			os.makedirs(outdir)

		for i in range(len(self.observationunit_list)):
			print('number of data:%d' % (i+1))
			event_list = []
			for obsunit in self.observationunit_list[0:i+1]:
				event_list.append(obsunit.param['nievt_cl_bary_phase_grpflag_esel'])
			print(event_list)
			xrayevtlist = nicrabgrp.xray.XrayEventList(event_list)
			outprofilefits = '%s/niacc%02d/niacc%02d_profile.fits' % (outdir,i+1,i+1)
			xrayevtlist.generate_profile_fitsfile(outprofilefits)
			for nphase in nicrabgrp.xray.NPHASE_LIST_DEFAULT:
					outpdf = '%s_n%02d.pdf' % (outprofilefits.replace('.fits',''),nphase)
					title  = 'Sum:%d, Phase N%d, ' % (i+1,nphase)
					title += '{:.1e} X-rays, {:.1e} MP-GRPs ({:.3f}%)'.format(xrayevtlist.numofevt_all,xrayevtlist.numofevt_mpgrp,xrayevtlist.fraction_ipgrp*100.0)
					xrayprofile = nicrabgrp.xray.XrayProfile(outprofilefits)

					enhancement, significance = xrayprofile.get_enhancement_significance(
						target_bins_list=self.param['GRP_SEARCH_TARGET_BINS_LIST'],
						nphase=self.param['GRP_SEARCH_NPHASE'])
					legend_title = 'Enhance:%.2f%% (%.2f-sigma)' % ((enhancement-1.0)*100.0,significance)
					xrayprofile.plot_profile_fitsfile(nphase=nphase,outpdf=outpdf,xmin=0.0,xmax=2.0,ymin=None,ymax=None,title=title,legend_title=legend_title)
					outpdf_zoom = '%s_n%02d_zoom.pdf' % (outprofilefits.replace('.fits',''),nphase)			
					xrayprofile.plot_profile_fitsfile(nphase=nphase,outpdf=outpdf_zoom,xmin=0.94,xmax=1.06,ymin=None,ymax=None,title=title,legend_title=legend_title)

