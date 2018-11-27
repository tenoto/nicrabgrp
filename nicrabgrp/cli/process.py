#!/usr/bin/env python
"""
HISTORY
2018-11-27 transfered from giantradiopulse to nicrabgrp library
2018-10-24 created by T.Enoto 
"""

import click 
import nicrabgrp.process

@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx):
	if ctx.invoked_subcommand is None:
		print(ctx.get_help())

@cli.command(help="Show parameters.")
@click.argument("file_path",type=click.Path(exists=True))
def show_parameters(file_path):
	process_manager = nicrabgrp.process.ProcessManager(file_path)

@cli.command(help="Convert radio text files to fitsfiles.")
@click.argument("file_path", type=click.Path(exists=True))
@click.option("--outdir", type=click.Path(), default='out/crabgrp')
def convert_radiofiles(file_path,outdir):
	process_manager = nicrabgrp.process.ProcessManager(file_path,outdir=outdir)
	process_manager.convert_radiofiles()

@cli.command(help="Prepare x-ray files: barycentric correction and adding PULSE_PHASE.")
@click.argument("file_path", type=click.Path(exists=True))
@click.option("--outdir", type=click.Path(), default='out/crabgrp')
def prepare_xrayfiles(file_path,outdir):
	process_manager = nicrabgrp.process.ProcessManager(file_path,outdir=outdir)
	process_manager.prepare_xrayfiles()

@cli.command(help="Adding GRP flags to event fits file.")
@click.argument("file_path", type=click.Path(exists=True))
@click.option("--outdir", type=click.Path(), default='out/crabgrp')
def add_grpflag_to_xrayfiles(file_path,outdir):
	process_manager = nicrabgrp.process.ProcessManager(file_path,outdir=outdir)
	process_manager.add_grpflag_to_xrayfiles()	

def main():
	cli()

if __name__ == "__main__":
	main()
