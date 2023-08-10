#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 14:15:14 2023

@author: kacharov
"""
import pytest
from astropy import units as u
from hawki_etc import hawki_etc


@pytest.fixture
def etc():
    return hawki_etc()


def test_snr2texp(etc):
    # Test the snr2texp method
    texp1 = etc.snr2texp(mag=20, band='Ks', snr=10,
                         airmass=1, seeing=0.8*u.arcsec, dit=10*u.s)
    texp2 = etc.snr2texp(mag=19*u.mag, band='Ks', snr=10,
                         airmass=1, seeing=0.8*u.arcsec, dit=10*u.s)
    assert (texp1 > 0 * u.s) & (
        texp2 > 0 * u.s), "snr2texp should return a positive exposure time"
    assert (texp1 > texp2), "The fainter star should return longer exposure time"

    texp1 = etc.snr2texp(mag=20, band='Ks', snr=15,
                         airmass=1, seeing=0.8*u.arcsec, dit=10*u.s)
    texp2 = etc.snr2texp(mag=20*u.mag, band='Ks', snr=10,
                         airmass=1, seeing=0.8*u.arcsec, dit=10*u.s)
    assert (
        texp1 > texp2), "The star with higher SNR requirement should return longer exposure time"

    texp1 = etc.snr2texp(mag=20, band='Ks', snr=10,
                         airmass=1.5, seeing=0.8*u.arcsec, dit=10*u.s)
    texp2 = etc.snr2texp(mag=20*u.mag, band='Ks', snr=10,
                         airmass=1, seeing=0.8*u.arcsec, dit=10*u.s)
    assert (
        texp1 > texp2), "The star with lower airmass should return longer exposure time"

    texp1 = etc.snr2texp(mag=20, band='Ks', snr=10,
                         airmass=1.5, seeing=1.0*u.arcsec, dit=10*u.s)
    texp2 = etc.snr2texp(mag=20*u.mag, band='Ks', snr=10,
                         airmass=1.5, seeing=0.8*u.arcsec, dit=10*u.s)
    assert (
        texp1 > texp2), "The star with worse seeing should return longer exposure time"

    texp1 = etc.snr2texp(mag=20, band='Ks', snr=10,
                         airmass=1.5, seeing=0.8*u.arcsec, dit=10*u.s)
    texp2 = etc.snr2texp(mag=20*u.mag, band='Ks', snr=10,
                         airmass=1.5, seeing=0.8*u.arcsec, dit=15*u.s)
    assert (
        texp1 > texp2), "The star with shorter dit should return longer exposure time"


def test_texp2snr(etc):
    # Test the texp2snr method
    snr1 = etc.texp2snr(mag=19, band='Ks', texp=100*u.s,
                        airmass=1.2, seeing=0.8*u.arcsec, dit=10*u.s)
    snr2 = etc.texp2snr(mag=20*u.mag, band='Ks', texp=100 *
                        u.s, airmass=1.2, seeing=0.8*u.arcsec, dit=10*u.s)
    assert (snr1 > 0) & (snr2 > 0), "texp2snr should return a positive SNR"
    assert (snr1 > snr2), "The fainter star should return lower SNR"

    snr1 = etc.texp2snr(mag=20, band='Ks', texp=150*u.s,
                        airmass=1.2, seeing=0.8*u.arcsec, dit=10*u.s)
    snr2 = etc.texp2snr(mag=20*u.mag, band='Ks', texp=100 *
                        u.s, airmass=1.2, seeing=0.8*u.arcsec, dit=10*u.s)
    assert (snr1 > snr2), "The star with shorter exposure should return lower SNR"

    snr1 = etc.texp2snr(mag=20, band='Ks', texp=100*u.s,
                        airmass=1.2, seeing=0.8*u.arcsec, dit=10*u.s)
    snr2 = etc.texp2snr(mag=20*u.mag, band='Ks', texp=100 *
                        u.s, airmass=1.8, seeing=0.8*u.arcsec, dit=10*u.s)
    assert (snr1 > snr2), "The star with higher airmass should return lower SNR"

    snr1 = etc.texp2snr(mag=20, band='Ks', texp=100*u.s,
                        airmass=1.2, seeing=0.8*u.arcsec, dit=10*u.s)
    snr2 = etc.texp2snr(mag=20*u.mag, band='Ks', texp=100 *
                        u.s, airmass=1.2, seeing=1.2*u.arcsec, dit=10*u.s)
    assert (snr1 > snr2), "The star with worse seeing should return lower SNR"

    snr1 = etc.texp2snr(mag=20, band='Ks', texp=100*u.s,
                        airmass=1.2, seeing=0.8*u.arcsec, dit=20*u.s)
    snr2 = etc.texp2snr(mag=20*u.mag, band='Ks', texp=100 *
                        u.s, airmass=1.2, seeing=0.8*u.arcsec, dit=10*u.s)
    assert (snr1 > snr2), "The star shorter DIT should return lower SNR"
