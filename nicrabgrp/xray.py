__author__  = 'Teruaki Enoto'
__version__ = '0.03'
__date__    = '2018 November 28'
"""
HISTORY
2018-11-27 mofified from xrayevent.py  
2018-11-15 created by T.Enoto 
"""

import os 
import sys
import numpy as np 
import astropy.io.fits as fits

import matplotlib.pyplot as plt 
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.size'] = '14'
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['xtick.top'] = 'True'
mpl.rcParams['ytick.right'] = 'True'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#mpl.rcParams['axes.grid'] = 'True'
mpl.rcParams['axes.xmargin'] = '.05' #'.05'
mpl.rcParams['axes.ymargin'] = '.05'
mpl.rcParams['savefig.facecolor'] = 'None'
mpl.rcParams['savefig.edgecolor'] = 'None'
mpl.rcParams['savefig.bbox'] = 'tight'

NPHASE_DEFAULT = 250 
NPHASE_LIST_DEFAULT = [50,100,200,250]
COLNAME_MPGRPFLAG_DEFAULT = 'MPGRP_SN5.0'
COLNAME_IPGRPFLAG_DEFAULT = 'IPGRP_SN5.0'

def get_profile(phase_array,nphase=NPHASE_DEFAULT):
	nevt = len(phase_array)
	profile_count, profile_binedge, profile_patches = plt.hist(phase_array,nphase,range=[0.0,1.0],histtype='step')
	profile_bincenter = 0.5*(profile_binedge[1:]+profile_binedge[:-1])	
	x = np.array(profile_bincenter)
	y = np.array(profile_count)
	xe = np.full(len(x),0.5/float(nphase))
	ye = np.array(np.sqrt(profile_count))
	ny = y / float(nevt)
	nye = ye / float(nevt)
	return x,xe,y,ye,ny,nye

def plot_two_profiles(x1,x1e,y1,y1e,x2,x2e,y2,y2e,outpdf,
		label1='',label2='',ylabel='',title='',legend_title='',
		xmin=0.0,xmax=2.0,ymin=None,ymax=None):

	outdir = os.path.dirname(outpdf)
	if outdir != '' and not os.path.exists(outdir):
		os.makedirs(outdir)

	plt.clf()
	fig, axes = plt.subplots(1,1,figsize=(9.6,5.6))
	plt.errorbar(np.hstack((x1,x1+1.0)),np.tile(y1,2),xerr=np.tile(x1e,2),yerr=np.tile(y1e,2),
		marker='',color='k',drawstyle='steps-mid',linewidth=1.0,label=label1)
	plt.errorbar(np.hstack((x2,x2+1.0)),np.tile(y2,2),xerr=np.tile(x2e,2),yerr=np.tile(y2e,2),
		marker='',color='r',drawstyle='steps-mid',linewidth=1.0,label=label2)	
	if ymax != None:
		plt.vlines([1.0],0.0,ymax,'k',linestyles='dashed',linewidth=1.0)  
	plt.title(title)
	legend = axes.legend(loc='upper right',shadow=False,fontsize=11.0,title=legend_title)
	axes.set_xlim(float(xmin),float(xmax))
	if ymin != None and ymax != None:
		axes.set_ylim(float(ymin),float(ymax))
	axes.set_xlabel('Pulse Phase')			
	axes.set_ylabel(ylabel)
	axes.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.savefig(outpdf,dpi=300)

def generate_profile_fitsfile(all_pulse_phase_array,mpgrp_flag_array,ipgrp_flag_array,outfitsfile,nphase_list=NPHASE_LIST_DEFAULT):
	print('--- method: %s ---' % sys._getframe().f_code.co_name)

	outdir = os.path.dirname(outfitsfile)
	if outdir != '' and not os.path.exists(outdir):
		os.makedirs(outdir)

	hdu_list = []
	hdu_list.append(fits.PrimaryHDU())

	cols_list = []
	for nphase in nphase_list:
		cols_list.append([])
		all_phase,all_phase_error,all_count,all_count_error,all_norm,all_norm_error = get_profile(all_pulse_phase_array,nphase)	
		cols_list[-1].append(fits.Column(name='PULSE_PHASE',format='1D',array=all_phase))		
		cols_list[-1].append(fits.Column(name='PHASE_ERROR',format='1D',array=all_phase_error))		

		cols_list[-1].append(fits.Column(name='ALL_COUNTS',format='K',array=all_count,unit='counts'))
		cols_list[-1].append(fits.Column(name='ALL_COUNTS_ERROR',format='1D',array=all_count_error,unit='counts'))	
		cols_list[-1].append(fits.Column(name='ALL_NORM_COUNTS',format='1D',array=all_norm,unit='au'))
		cols_list[-1].append(fits.Column(name='ALL_NORM_COUNTS_ERROR',format='1D',array=all_norm_error,unit='au'))		

		mpgrp_phase,mpgrp_phase_error,mpgrp_count,mpgrp_count_error,mpgrp_norm,mpgrp_norm_error = get_profile(all_pulse_phase_array[mpgrp_flag_array],nphase)	
		cols_list[-1].append(fits.Column(name='MPGRP_COUNTS',format='K',array=mpgrp_count,unit='counts'))
		cols_list[-1].append(fits.Column(name='MPGRP_COUNTS_ERROR',format='1D',array=mpgrp_count_error,unit='counts'))	
		cols_list[-1].append(fits.Column(name='MPGRP_NORM_COUNTS',format='1D',array=mpgrp_norm,unit='au'))
		cols_list[-1].append(fits.Column(name='MPGRP_NORM_COUNTS_ERROR',format='1D',array=mpgrp_norm_error,unit='au'))		

		nonmpgrp_phase,nonmpgrp_phase_error,nonmpgrp_count,nonmpgrp_count_error,nonmpgrp_norm,nonmpgrp_norm_error = get_profile(all_pulse_phase_array[~mpgrp_flag_array],nphase)	
		cols_list[-1].append(fits.Column(name='NONMPGRP_COUNTS',format='K',array=nonmpgrp_count,unit='counts'))
		cols_list[-1].append(fits.Column(name='NONMPGRP_COUNTS_ERROR',format='1D',array=nonmpgrp_count_error,unit='counts'))	
		cols_list[-1].append(fits.Column(name='NONMPGRP_NORM_COUNTS',format='1D',array=nonmpgrp_norm,unit='au'))
		cols_list[-1].append(fits.Column(name='NONMPGRP_NORM_COUNTS_ERROR',format='1D',array=nonmpgrp_norm_error,unit='au'))		

		ipgrp_phase,ipgrp_phase_error,ipgrp_count,ipgrp_count_error,ipgrp_norm,ipgrp_norm_error = get_profile(all_pulse_phase_array[ipgrp_flag_array],nphase)	
		cols_list[-1].append(fits.Column(name='IPGRP_COUNTS',format='K',array=ipgrp_count,unit='counts'))
		cols_list[-1].append(fits.Column(name='IPGRP_COUNTS_ERROR',format='1D',array=ipgrp_count_error,unit='counts'))	
		cols_list[-1].append(fits.Column(name='IPGRP_NORM_COUNTS',format='1D',array=ipgrp_norm,unit='au'))
		cols_list[-1].append(fits.Column(name='IPGRP_NORM_COUNTS_ERROR',format='1D',array=ipgrp_norm_error,unit='au'))		

		nonipgrp_phase,nonipgrp_phase_error,nonipgrp_count,nonipgrp_count_error,nonipgrp_norm,nonipgrp_norm_error = get_profile(all_pulse_phase_array[~ipgrp_flag_array],nphase)	
		cols_list[-1].append(fits.Column(name='NONIPGRP_COUNTS',format='K',array=nonipgrp_count,unit='counts'))
		cols_list[-1].append(fits.Column(name='NONIPGRP_COUNTS_ERROR',format='1D',array=nonipgrp_count_error,unit='counts'))	
		cols_list[-1].append(fits.Column(name='NONIPGRP_NORM_COUNTS',format='1D',array=nonipgrp_norm,unit='au'))
		cols_list[-1].append(fits.Column(name='NONIPGRP_NORM_COUNTS_ERROR',format='1D',array=nonipgrp_norm_error,unit='au'))		

		mpgrp_norm_sub = mpgrp_norm - nonmpgrp_norm
		mpgrp_norm_sub_error = np.sqrt(mpgrp_norm_error**2 + nonmpgrp_norm_error**2)
		cols_list[-1].append(fits.Column(name='MPGRP_NORM_SUB',format='1D',array=mpgrp_norm_sub,unit='au'))
		cols_list[-1].append(fits.Column(name='MPGRP_NORM_SUB_ERROR',format='1D',array=mpgrp_norm_sub_error,unit='au'))	

		extname = 'PROFILE_N%d' % nphase
		hdu_list.append(fits.BinTableHDU.from_columns(cols_list[-1],name=extname))

		hdu_list[-1].header['NPHASE'] = nphase 
		hdu_list[-1].header['NXALL']   = len(all_pulse_phase_array)
		hdu_list[-1].header['NXMPGRP'] = len(all_pulse_phase_array[mpgrp_flag_array])		
		hdu_list[-1].header['NXIPGRP'] = len(all_pulse_phase_array[ipgrp_flag_array])				
		hdu_list[-1].header['NXNMPGRP'] = len(all_pulse_phase_array[~mpgrp_flag_array])		
		hdu_list[-1].header['NXNIPGRP'] = len(all_pulse_phase_array[~ipgrp_flag_array])		
		hdu_list[-1].header['FR_MPGRP'] = float(len(all_pulse_phase_array[mpgrp_flag_array]))/float(len(all_pulse_phase_array))
		hdu_list[-1].header['FR_IPGRP'] = float(len(all_pulse_phase_array[ipgrp_flag_array]))/float(len(all_pulse_phase_array))
	fits.HDUList(hdu_list).writeto(outfitsfile,overwrite=True)	


class XrayEventFile():
	def __init__(self,eventfits):
		self.eventfits = eventfits 

		if not os.path.exists(self.eventfits):
			sys.stderr.write('file %s does not exist.\n' % self.eventfits)
			quit()
		print('XrayEvent: {} loaded.'.format(self.eventfits))
		print('XrayEvent: {} GB'.format(1e-3*(os.path.getsize(self.eventfits)>>20)))

		try:
			with fits.open(self.eventfits) as self.hdu:
				self.header  = self.hdu['EVENTS'].header				
				self.data    = self.hdu['EVENTS'].data
				self.columns = self.hdu['EVENTS'].columns
		except OSError as e:
			raise 
		self.numofevt_all = len(self.data['TIME'])
		print('XrayEvent: accessed to event list ({:,} events)'.format(self.numofevt_all))

	def show_statistics(self):
		print('numofevt_all: {:,} events'.format(self.numofevt_all))

	def add_grpflag(self,mpgrplist_fitsfile,ipgrplist_fitsfile,grp_snthr,outfitsfile,out_column_name_mpgrp,out_column_name_ipgrp):
		print('--- %s method: %s ---' % (self.eventfits,sys._getframe().f_code.co_name))		
		print("MP-GRP list fitsfile: {}".format(mpgrplist_fitsfile))
		print("IP-GRP list fitsfile: {}".format(ipgrplist_fitsfile))		
		print("GRP SN threshold: {}".format(grp_snthr))		
		print("output fitsfile: {}".format(outfitsfile))		
		print("output column name MP-GRP: {}".format(out_column_name_mpgrp))				
		print("output column name IP-GRP: {}".format(out_column_name_ipgrp))	

		### Radio GRP MP
		if not os.path.exists(mpgrplist_fitsfile):
			sys.stderr.write('file %s does not exist.\n' % mpgrplist_fitsfile)
			quit()
		try:
			with fits.open(mpgrplist_fitsfile) as hdu_mpgrp:
				hdu_mpgrp_table = hdu_mpgrp['GRP'].data
		except OSError as e:
			raise 

		flag_mpgrp_snthr_filter = hdu_mpgrp_table['PEAK_SN'] > grp_snthr
		mpgrp_mod_pulse_number_filtered = hdu_mpgrp_table['MOD_PULSE_NUMBER'][flag_mpgrp_snthr_filter]

		flag_xray_isin_mpgrp = np.isin(self.data['MOD_PULSE_NUMBER'],mpgrp_mod_pulse_number_filtered)
		new_column_flag_xray_isin_mpgrp = np.array(flag_xray_isin_mpgrp,dtype='bool')	

		print(hdu_mpgrp_table['MOD_PULSE_NUMBER'],len(hdu_mpgrp_table['MOD_PULSE_NUMBER']))
		print(mpgrp_mod_pulse_number_filtered,len(mpgrp_mod_pulse_number_filtered))
		print(flag_xray_isin_mpgrp,len(flag_xray_isin_mpgrp))

		### Radio GRP IP
		if not os.path.exists(ipgrplist_fitsfile):
			sys.stderr.write('file %s does not exist.\n' % ipgrplist_fitsfile)
			quit()
		try:
			with fits.open(ipgrplist_fitsfile) as hdu_ipgrp:
				hdu_ipgrp_table = hdu_ipgrp['GRP'].data
		except OSError as e:
			raise 
			
		flag_ipgrp_snthr_filter = hdu_ipgrp_table['PEAK_SN'] > grp_snthr
		ipgrp_mod_pulse_number_filtered = hdu_ipgrp_table['MOD_PULSE_NUMBER'][flag_ipgrp_snthr_filter]

		flag_xray_isin_ipgrp = np.isin(self.data['MOD_PULSE_NUMBER'],ipgrp_mod_pulse_number_filtered)
		new_column_flag_xray_isin_ipgrp = np.array(flag_xray_isin_ipgrp,dtype='bool')	
		new_columns = fits.ColDefs([
			fits.Column(name=out_column_name_mpgrp,format='L',array=new_column_flag_xray_isin_mpgrp),
			fits.Column(name=out_column_name_ipgrp,format='L',array=new_column_flag_xray_isin_ipgrp)			
			])
			
		hdu_primary = fits.PrimaryHDU()
		hdu_newtable = fits.BinTableHDU.from_columns(self.columns+new_columns,name='EVENTS')	
		hdulist = fits.HDUList([hdu_primary,hdu_newtable])
		hdulist.writeto(outfitsfile,overwrite=True)	

		cmd = 'ftappend %s+2 %s' % (self.eventfits,outfitsfile)
		print(cmd);os.system(cmd)

		cmd = 'cphead %s+0 %s+0;\n' % (self.eventfits,outfitsfile)
		cmd += 'cphead %s+1 %s+1;\n' % (self.eventfits,outfitsfile)		
		print(cmd);os.system(cmd)

class XrayEventList():
	def __init__(self,evtfile_list):
		self.evtfile_list = evtfile_list
		print('evtfile_list: {}'.format(self.evtfile_list))
		self.set_xrayevt_list()
		self.show_statistics()

	def set_xrayevt_list(self,colname_mpgrpflag=COLNAME_MPGRPFLAG_DEFAULT,colname_ipgrpflag=COLNAME_IPGRPFLAG_DEFAULT):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)		
		self.hdu_list = []
		for evtfile in self.evtfile_list:
			self.hdu_list.append(XrayEventFile(evtfile))

		self.time_array_list = []
		self.pulse_phase_array_list = []
		self.mpgrp_flag_array_list = []
		self.ipgrp_flag_array_list = []
		for hdu in self.hdu_list:
			self.time_array_list.append(hdu.data['TIME'])
			self.pulse_phase_array_list.append(hdu.data['PULSE_PHASE'])
			self.mpgrp_flag_array_list.append(hdu.data[COLNAME_MPGRPFLAG_DEFAULT])			
			self.ipgrp_flag_array_list.append(hdu.data[COLNAME_IPGRPFLAG_DEFAULT])					

		self.time_array = np.hstack(self.time_array_list)
		self.pulse_phase_array = np.hstack(self.pulse_phase_array_list) 
		self.mpgrp_flag_array = np.hstack(self.mpgrp_flag_array_list)
		self.ipgrp_flag_array = np.hstack(self.ipgrp_flag_array_list)

		self.numofevt_all = len(self.time_array)
		self.numofevt_mpgrp = len(self.mpgrp_flag_array[self.mpgrp_flag_array])
		self.numofevt_ipgrp = len(self.ipgrp_flag_array[self.ipgrp_flag_array])	
		self.fraction_mpgrp = float(self.numofevt_mpgrp)/float(self.numofevt_all)
		self.fraction_ipgrp = float(self.numofevt_ipgrp)/float(self.numofevt_all)	

	def show_statistics(self):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)	
		print('numofevt_all: {:,}'.format(self.numofevt_all))
		print('numofevt_mpgrp: {:,} ({:.2f}%)'.format(self.numofevt_mpgrp,self.fraction_mpgrp*100.0))		
		print('numofevt_ipgrp: {:,} ({:.2f}%)'.format(self.numofevt_ipgrp,self.fraction_ipgrp*100.0))				

	def generate_profile_fitsfile(self,outfitsfile,nphase_list=NPHASE_LIST_DEFAULT):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)

		generate_profile_fitsfile(
			all_pulse_phase_array=self.pulse_phase_array,
			mpgrp_flag_array=self.mpgrp_flag_array,
			ipgrp_flag_array=self.ipgrp_flag_array,
			outfitsfile=outfitsfile,
			nphase_list=nphase_list)

class XrayProfile():
	def __init__(self,fitsfile):
		self.fitsfile = fitsfile 

		if not os.path.exists(self.fitsfile):
			sys.stderr.write('file %s does not exist.\n' % self.fitsfile)
			quit()
		self.hdu = fits.open(self.fitsfile)
		print('XrayProfile: {} loaded.'.format(self.fitsfile))

	def plot_profile_fitsfile(self,nphase,outpdf,xmin=0.0,xmax=2.0,ymin=None,ymax=None,title='',legend_title=''):
		print('--- %s method: %s ---' % (self.fitsfile,sys._getframe().f_code.co_name))

		extname = 'PROFILE_N%d' % nphase 
		x1  = self.hdu[extname].data['PULSE_PHASE']
		x1e = self.hdu[extname].data['PHASE_ERROR']
		y1  = self.hdu[extname].data['ALL_NORM_COUNTS']
		y1e = self.hdu[extname].data['ALL_NORM_COUNTS_ERROR']
		y2  = self.hdu[extname].data['MPGRP_NORM_COUNTS']
		y2e = self.hdu[extname].data['MPGRP_NORM_COUNTS_ERROR']

		plot_two_profiles(x1,x1e,y1,y1e,x1,x1e,y2,y2e,outpdf,
			label1='All',label2='MPGRP',ylabel='Normalized counts',title=title,legend_title=legend_title,
			xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)

	def get_enhancement_significance(self,target_bins_list,nphase):
		extname = 'PROFILE_N%d' % nphase 
		data = self.hdu[extname].data

		nongrp_sum_counts = sum(data['NONMPGRP_COUNTS'])
		grp_sum_counts = sum(data['MPGRP_COUNTS'])

		nongrp_bottom_counts = min(data['NONMPGRP_COUNTS'])
		nongrp_bottom_norm = float(nongrp_bottom_counts) / float(nongrp_sum_counts)

		grp_enhance_counts = 0
		nongrp_enhance_counts = 0			
		for target_bin in target_bins_list:
			grp_enhance_counts += data['MPGRP_COUNTS'][target_bin]
			nongrp_enhance_counts += data['NONMPGRP_COUNTS'][target_bin]
		grp_enhance_norm = float(grp_enhance_counts) / float(grp_sum_counts)
		nongrp_enhance_norm = float(nongrp_enhance_counts) / float(nongrp_sum_counts)

		grp_enhance_error_norm = np.sqrt(float(grp_enhance_counts)) / float(grp_sum_counts)
		nongrp_enhance_error_norm = np.sqrt(float(nongrp_enhance_counts)) / float(nongrp_sum_counts)

		enhancement  = (grp_enhance_norm - nongrp_bottom_norm) / (nongrp_enhance_norm - nongrp_bottom_norm)
		significance = (grp_enhance_norm - nongrp_enhance_norm) / np.sqrt(grp_enhance_error_norm**2+nongrp_enhance_error_norm**2)
		return enhancement,significance


