#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 11:58:54 2023

@author: kacharov
"""

from hawki_etc import hawki_etc
from astropy import units as u

etc = hawki_etc()

texp = etc.snr2texp(mag=22.5*u.mag, band='Ks', snr=5, airmass=1.2,
                    seeing=0.8*u.arcsec, dit=10*u.s)

snr = etc.texp2snr(mag=22.5, band='Ks', texp=3600*u.s,
                   airmass=1.2, seeing=0.8*u.arcsec, dit=10*u.s)

print()
print('Texp = ', texp)
print('SNR = ', snr)
print()
print(dir(etc))