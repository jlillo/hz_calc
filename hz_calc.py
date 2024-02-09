import sys
import os.path

import numpy as np
from math import *
from astropy import constants as c
import astropy.units as u
import argparse

"""HZ waraper
	
	Wraper from the Ravi Kopparapu functions to obtain the HZ period/sma for a given set of 
	provided stellar and orbital values. 

	INPUT
	-----
	Teff		float		: Stellar effective temperature

	OPTIONAL
	--------
	--PER		boolean		: provide the result in days, default=False
	--MS		float		: stellar mass if period=True

	OUTPUT
	------
	opt_in, opt_out			: interior and exterior HZ limits for the OPTIMISTIC scenario
	cons_in, cons_out		: interior and exterior HZ limits for the CONSERVATIVE scenario
"""


def cli():
    """command line inputs

    Get parameters from command line

    Returns
    -------
    Arguments passed by command line
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("Teff", help="Effective temperature")
    parser.add_argument("-P", "--PER", help="Output in days (assuming MS)", action="store_true")
    parser.add_argument("-M", "--MS", help="Stellar mass")
    parser.add_argument("-V", "--VERBOSE", help="Print the results in terminal",action="store_true")
    args = parser.parse_args()
    return args


def hz(teff,Lum,DO_SEFF=False):

	seff = [0,0,0,0,0,0]
	seffsun  = [1.776,1.107, 0.356, 0.320, 1.188, 0.99]
	a = [2.136e-4, 1.332e-4, 6.171e-5, 5.547e-5, 1.433e-4, 1.209e-4]
	b = [2.533e-8, 1.580e-8, 1.698e-9, 1.526e-9, 1.707e-8, 1.404e-8]
	c = [-1.332e-11, -8.308e-12, -3.198e-12, -2.874e-12, -8.968e-12, -7.418e-12]
	d = [-3.097e-15, -1.931e-15, -5.575e-16, -5.011e-16, -2.084e-15, -1.713e-15]

	tstar = teff - 5780.
	for i in range(len(a)):
		seff[i] = seffsun[i] + a[i]*tstar + b[i]*tstar**2 + c[i]*tstar**3 + d[i]*tstar**4

	hz_au = np.sqrt(Lum/seff)

	result = hz_au
	if DO_SEFF:
		result = seff

	return result

def myLogFormat(y,pos):
    # Find the number of decimal places required
    decimalplaces = int(np.maximum(-np.log10(y),0))     # =0 for numbers >=1
    # Insert that number into a format string
    formatstring = '{{:.{:1d}f}}'.format(decimalplaces)
    # Return the formatted tick label
    return formatstring.format(y)

def P_from_sma(sma,Ms,VERBOSE=False):
	"""
		Estimate Planet period from semi-major axis and stellar mass
		using the 3rd Kepler law.
		Input: sma [au], Ms [Msun]
	"""
	period = np.sqrt(((sma*c.au.value)**3 *4*np.pi**2) / (c.G.value*Ms*c.M_sun.value) ) / (3600*24)
	if VERBOSE: print("Period [days] = ", period)
	return period


def get_hz(Teff,PER=False,Ms=None):

	teff = float(Teff)


	# =====  HABITABLE ZONE calculations

	tab = np.genfromtxt('isocz019_SolarAge.dat',dtype=None) #logage,miso,mact,logl,logtiso,loggiso
	logtiso = tab[:,4]
	logl = tab[:,3]
	tiso = 10.**logtiso[0:41]
	Lumiso = 10.**logl[0:41]
	Mact = tab[0:41,2]

	LL = np.interp(teff, tiso, Lumiso)
	res = hz(teff,LL, DO_SEFF=False)

	opt_in   = res[0] # recent Venus
	opt_out  = res[3] # early Mars
	cons_in  = res[4] # 5 Mearth Runaway Greenhouse.
	cons_out = res[2] # max Greenhouse

	if PER:
		if Ms == None:
			print("\t --> Estimating the Mstar with the isochrone isocz019_SolarAge.dat")
			Mstar = np.interp(teff, tiso, Mact)
			print("\t     - The estimated mass is: ",np.round(Mstar,2)," Msun \n")
		else:
			Mstar = float(Ms)
		opt_in   = P_from_sma(opt_in,Mstar)
		opt_out  = P_from_sma(opt_out,Mstar)
		cons_in  = P_from_sma(cons_in,Mstar)
		cons_out = P_from_sma(cons_out,Mstar)

	return opt_in, opt_out, cons_in, cons_out


# ======================================
# 	        MAIN
# ======================================

if __name__ == "__main__":

	print("................................")
	print("    HZ calculator wraper        ")
	print("................................")

	args = cli()

	opt_in, opt_out, cons_in, cons_out = get_hz(args.Teff,PER=args.PER,Ms=args.MS)

	units = 'au'
	if args.PER: units = 'days'

	if args.VERBOSE:
		print("\t ====> Results:")
		print("\t HZ conservative ("+units+") = ", round(opt_in,2), "-->", round(opt_out,2))
		print("\t HZ optimistic ("+units+") = ", round(cons_in,2), "-->", round(cons_out,2))
	else:
		print("Use --VERBOSE option to print out the results")
	
	print("\n")
	# return opt_in, opt_out, cons_in, cons_out