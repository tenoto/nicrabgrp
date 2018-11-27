__author__  = 'Teruaki Enoto'
__version__ = '0.02'
__date__    = '2018 November 15'
"""
HISTORY
2018-11-27 mofified from xrayevent.py  
2018-11-15 created by T.Enoto 
"""

import os 
import sys
import numpy as np 
import astropy.io.fits as fits

NPHASE_DEFAULT = 250 

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

class XrayEventFile():
	def __init__(self,eventfits):
		self.eventfits = eventfits 

		if not os.path.exists(self.eventfits):
			sys.stderr.write('file %s does not exist.\n' % self.eventfits)
			quit()
		print('XrayEvent: {} loaded.'.format(self.eventfits))

		try:
			with fits.open(self.eventfits) as self.hdu:
				self.header  = self.hdu['EVENTS'].header				
				self.data    = self.hdu['EVENTS'].data
				self.columns = self.hdu['EVENTS'].columns
		except OSError as e:
			raise 

		self.numofevt_all = len(self.data['TIME'])

	def show_statistics(self):
		print('numofevt_all: {:,}'.format(self.numofevt_all))

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



