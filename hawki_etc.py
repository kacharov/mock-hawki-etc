#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 15:12:19 2023

@author: kacharov
"""

import numpy as np
from astropy import units as u
import hmbp


class hawki_etc:
    '''
    This is the HAWK-I ETC class.

    Parameters
    ----------
    None

    Attributes
    ----------
    diameter : Quantity
        Main mirror diameter of the VLT
    fov : Quantity
        Size of the focal plane (HAWK-I User manual)
    ron : int
        Readout noise per pixel
    dark : Quantity
        Dark current noise per pixel
    pixelscale : Quantity
        Pixel scale per pixel
    efficiency : float
        Fraction of photons converted to CCD counts.
    ext_coef : float
        Extinction coefficient.
    '''

    def __init__(self):
        # HAWK-I basic parameters

        self.diameter = 8*u.m  # main mirror diameter of the VLT

        # size of the focal plane (HAWK-I User manual)
        self.fov = 10.5*u.arcmin

        self.ron = 5  # readout noise - counts per px
        self.dark = 0.01*u.s**(-1)  # dark current counts per px
        self.pixelscale = 0.1064*u.arcsec  # per px

        # fraction of photons converted to CCD counts - this variable
        # incorporates loss of photons due to multiple factors:
        # CCD response rate, losses and holes in the mirrors, etc.
        # It is adjusted to give a limiting magnitude of Ks=22.5
        # for 1h exposure and SNR~5, at airmass=1 and seeing=0.8".
        self.efficiency = 0.85

        # assume exponential decrease of the number of photons with airmass
        # with an extinction coef. 0.4 (this factor depends on the pass band,
        # which we currently ignore).
        self.ext_coef = 0.4

    def get_counts(self, mag=0, band='H', airmass=1):
        '''
        This method is used to calculate the number of counts that hit the
        detector per second.

        Parameters
        ----------
        mag : float or quantity.Quantity, optional
            The magnitude of the star. If it is a float, Vega magnitudes are
            assumed. Default value is 0.0.
        band : str, optional
            The filter band. Default value is 'H'.
        airmass : float, optional
            The airmass of the observation. Default value is 1.0.

        Return
        ------
        counts : float
            The number of counts that hit the detector per second.
        '''

        # mag needs to have a magnitude unit. If it is just a number,
        # assume Vega magnitudes.
        if not isinstance(mag, u.quantity.Quantity):
            mag *= u.mag

        # get the number of photons per s per m2 in vaccuum
        n_phot = hmbp.for_flux_in_filter(band, mag,
                                         instrument="HAWKI",
                                         observatory="Paranal")

        # get the number of photons that hit the telescope per s
        # assume exponential decrease of the number of photons with airmass
        n_phot *= np.exp(-self.ext_coef*airmass) * np.pi * (self.diameter/2)**2

        counts = n_phot * self.efficiency

        return counts

    def get_sky(self, band='H', airmass=1, seeing=0.8*u.arcsec):
        '''
        This method is used to get the sky counts in the area of the detector
        occupied by the targtet, given a seeing condition.

        Parameters
        ----------
        band : str, optional
            The filter band. Default value is 'H'.
        airmass : float, optional
            The airmass of the observation. Default value is 1.0.
        seeing : quantity.Quantity, optional
            The seeing of the observation, expressed in Angular units (arcsec).
            Default value is 0.8*arcsec.

        Returns
        -------
        sky_counts : float
            The sky counts that hit the area of the target star.
        '''

        # get the sky counts per s per m2 in the focal plane
        n_sky = hmbp.in_skycalc_background(band, airmass=airmass,
                                           instrument="HAWKI",
                                           observatory="Paranal")

        # sky counts that hit the detector
        # assuming here that the detector covers completely the focal plane,
        # which is generally not true, but irrelevant for our purpose
        sky_counts = n_sky * np.pi * (self.diameter/2)**2 * self.efficiency

        # get the sky counts that hit the area of the target star,
        # given a seeing value
        sky_counts *= seeing**2 / (self.fov.to("arcsec")/2)**2

        return sky_counts

    def get_noise(self, counts, band='H', airmass=1, seeing=0.8*u.arcsec):
        '''
        This method is used to calculate the square of the Poisson noise in
        the area of the target on the detector, occupied by the target.
        It includes Poisson noise from the target, from the sky, dark current
        and readout noise. The latter is output separately, because it does
        not depend on the exposure time.

        Parameters
        ----------
        counts : float
            The number of counts that hit the telescope per second from the
            target.
        band : str, optional
            The filter band. Default value is 'H'.
        airmass : float, optional
            The airmass of the observation. Default value is 1.
        seeing : quantity.Quantity, optional
            The seeing of the observation, expressed in Angular units (arcsec).
            Default value is 0.8*arcsec.

        Returns
        -------
        noise_counts : float
            The Poisson noise per second.
        ron_counts : float
            The readout noise per DIT (Detector Integration Time).
        '''

        # get the sky counts per s
        self.sky_counts = self.get_sky(
            band=band, airmass=airmass, seeing=seeing)

        # compute the dark current noise in the area of the target
        dark_counts = self.dark * np.pi * seeing**2 / self.pixelscale**2 * u.ph

        # Poisson noise per s
        noise_counts = counts + 2*self.sky_counts + dark_counts

        # Poisson noise per DIT
        ron_counts = self.ron * np.pi * seeing**2 / self.pixelscale**2 * u.ph

        return noise_counts, ron_counts

    def snr2texp(self, mag=0, band='H', snr=5, airmass=1, seeing=0.8*u.arcsec,
                 dit=10*u.s):
        '''
        This method is used to calculate the exposure time for a given SNR
        (Signal-to-Noise Ratio).

        Parameters
        ----------
        mag : float or quantity.Quantity, optional
            The magnitude of the star. Default value is 0.
        band : str, optional
            The filter band. Default value is 'H'.
        snr : int, optional
            The signal-to-noise ratio. Default value is 5.
        airmass : float, optional
            The airmass of the observation. Default value is 1.
        seeing : quantity.Quantity, optional
            The seeing of the observation, expressed in Angular units (arcsec).
            Default value is 0.8*arcsec.
        dit : quantity.Quantity, optional
            The exposure time per frame. Default value is 10*seconds.

        Returns
        -------
        texp : float
            The total exposure time in seconds.
        '''

        # add the the user defined parametrs to the etc class
        self.mag = mag
        self.band = band
        self.snr = snr
        self.airmass = airmass
        self.seeing = seeing
        self.dit = dit

        # get the counts hitting the detector per s
        self.counts = self.get_counts(mag=mag, band=band, airmass=airmass)
        # estimate the Poisson noise per s
        noise, ron = self.get_noise(
            self.counts, band=band, airmass=airmass, seeing=seeing)

        # compute exposure time in seconds
        self.texp = (dit * snr**2 * noise + snr**2 * ron) / \
            (self.counts**2 * dit)
        self.texp *= u.ph

        self.ndit = np.round(self.texp / dit)

        return self.texp

    def texp2snr(self, mag=0, band='H', texp=10*u.s, airmass=1,
                 seeing=0.8*u.arcsec, dit=10*u.s):
        '''
        This method is used to calculate the signal-to-noise ratio for a given
        exposure time.

        Parameters
        ----------
        mag : float or quantity.Quantity, optional
            The magnitude of the star. Default value is 0.
        band : str, optional
            The filter band. Default value is 'H'.
        texp : quantity.Quantity, optional
            The total exposure time. Default value is 10*seconds.
        airmass : float, optional
            The airmass of the observation. Default value is 1.
        seeing : quantity.Quantity, optional
            The seeing of the observation, expressed in Angular units (arcsec).
            Default value is 0.8*arcsec.
        dit : quantity.Quantity, optional
            The exposure time per frame. Default value is 10*seconds.

        Returns
        -------
        snr : float
            The signal-to-noise ratio.
        '''

        # add the the user defined parametrs to the etc class
        self.mag = mag
        self.band = band
        self.texp = texp
        self.airmass = airmass
        self.seeing = seeing
        self.dit = dit
        self.ndit = np.round(texp / dit)

        # get the counts hitting the detector per s
        self.counts = self.get_counts(mag=mag, band=band, airmass=airmass)
        # estimate the Poisson noise per s
        noise, ron = self.get_noise(
            self.counts, band=band, airmass=airmass, seeing=seeing)

        self.snr = self.counts * texp / np.sqrt(noise * texp + self.ndit * ron)
        self.snr = self.snr.value

        return self.snr
