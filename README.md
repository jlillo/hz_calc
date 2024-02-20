# Habitable Zone calculator

This module is a wraper of the `hz.py` code from Ravi Koppparapu to estimate the optimistic and conservative habitable zone location for main-sequence stars.

### Usage

To get the HZ in AU units for a star with Teff=4500K:

`python hz_calc.py 4500`

If the user prefers to get the results in days:

`python hz_calc.py 4500 --VERBOSE --PER`

In this case, the code will estimate the mass of the star based on the solar metallicity isochrone. If you want to provide your stellar mass, then:

`python hz_calc.py 4500 --VERBOSE --PER --MS 0.75`

Importing it as a module

`import hz_calc`
`opt_in, opt_out, cons_in, cons_out = hz_calc.get_hz(4500,PER=true,Ms=0.75)`

### Citations needed

If you make use of this code, please do not forget to acknowledge the following references:
* Kopparapu et al., 2014, ApJ, 787, 29, *Habitable Zones around Main-sequence Stars: Dependence on Planetary Mass* [ADS link](https://ui.adsabs.harvard.edu/abs/2014ApJ...787L..29K/abstract)
* Lillo-Box et al., 2022, A&A, 667, 102, *The KOBE experiment: K-dwarfs Orbited By habitable Exoplanets. Project goals, target selection, and stellar characterization* [ADS link](https://ui.adsabs.harvard.edu/abs/2022A%26A...667A.102L/abstract)

