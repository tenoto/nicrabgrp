#!/usr/bin/env python
"""
HISTORY
2018-11-28 modified by T.Enoto
2018-10-09 modified by T.Enoto
2018-10-01 created by T.Enoto 
"""

import click 
import nicrabgrp.xray

@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx):
	if ctx.invoked_subcommand is None:
		print(ctx.get_help())

@cli.command(help="Generate X-ray pulse profile fitsfile.")
@click.argument("evtfiles", type=click.Path(exists=True),nargs=-1)
@click.option("--outfitsfile", type=click.Path(), default=None)
def generate_profile_fitsfile(evtfiles,outfitsfile):
	evtfile_list = []
	for evtfile in evtfiles:
		evtfile_list.append(evtfile)
	xrayevtlist = nicrabgrp.xray.XrayEventList(evtfile_list)
	xrayevtlist.generate_profile_fitsfile(outfitsfile)

@cli.command(help="Plot X-ray pulse profile fitsfile.")
@click.argument("fitsfile", type=click.Path(exists=True))
@click.option("--nphase", default=250)
@click.option("--outpdf",default='out.pdf')
@click.option("--xmin",default=0.0)
@click.option("--xmax",default=2.0)
@click.option("--ymin",default=None)
@click.option("--ymax",default=None)
def plot_profile_fitsfile(fitsfile,nphase,outpdf,xmin,xmax,ymin,ymax):
	xrayprofile = nicrabgrp.xray.XrayProfile(fitsfile)
	xrayprofile.plot_profile_fitsfile(nphase,outpdf,xmin,xmax,ymin,ymax)

@cli.command(help="Calculate enahancement significance. (target bin starts from 1 to 250)")
@click.argument("fitsfile", type=click.Path(exists=True))
@click.argument("target_bins")
@click.option("--nphase", default=250)
def get_enhancement_significance(fitsfile,target_bins,nphase):
	target_bins_list = []
	for target_bin_txt in target_bins.split(','):
		target_bins_list.append(int(target_bin_txt))
	xrayprofile = nicrabgrp.xray.XrayProfile(fitsfile)
	enhancement, significance = xrayprofile.get_enhancement_significance(target_bins_list,nphase)
	print('enhancement: {:.4f}'.format(enhancement))
	print('significance: {:.4f}'.format(significance))

def main():
	cli()

if __name__ == "__main__":
	main()
