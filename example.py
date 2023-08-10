#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 11:58:54 2023

@author: kacharov
"""

from hawki_etc import hawki_etc
from astropy import units as u

# Initiate the HAWK-I ETC
etc = hawki_etc()

# Compute the total exposure time given a set of parameters and a target SNR
texp = etc.snr2texp(mag=22.5*u.mag, band='Ks', snr=5, airmass=1.2,
                    seeing=0.8*u.arcsec, dit=10*u.s)

# Compute the SNR given a set of parameters and a total exposure time
snr = etc.texp2snr(mag=22.5, band='Ks', texp=3600*u.s,
                   airmass=1.2, seeing=0.8*u.arcsec, dit=10*u.s)

# Print the results
print()
print('Texp = ', texp)
print('SNR = ', snr)
print()

# Reveal all attributes of the etc class
print(dir(etc))