# This is Python 3 code but with these import there's at least a chance it'll work in Python 2
from __future__ import division, print_function, unicode_literals, absolute_import, with_statement

import numpy as np
import astropy.io.fits as fits
import astropy.units as u

class ZodiacalLight:
    """

    """
    # Colina, Bohlin & Castelli solar spectrum normalised to V band flux 
    # of 184.2 ergs/s/cm^2/A,
    solar_normalisation = 184.2 * u.erg * u.second**-1 * u.cm**-2 * u.Angstrom**-1
    # Leinert at al NEP zodical light 1.81e-18 erg/s/cm^2/A/arcsec^2 at 0.5 um, 
    # Aldering suggests using 0.01 dex lower.
    zl_nep = 1.81e-18 * u.erg * u.second**-1 * u.cm**-2 * u.Angstrom**-1 * \
             u.arcsecond**-2 * 10**(-0.01)
    zl_normalisation = zl_nep / solar_normalisation
    # Central wavelength
    lambda_c = 0.5 * u.micron
    # Aldering reddening parameters
    f_blue = 0.9
    f_red = 0.48    

    def __init__(self, solar_path='../resources/sun_castelli.fits'):
        # Pre-calculate zodiacal light spectrum for later use.
        self._calculate_spectrum(solar_path)

    def _calculate_spectrum(self, solar_path):
        """
        
        """
        # Load absolute solar spectrum from Collina, Bohlin & Castelli (1996)
        sun = fits.open(solar_path)
        sun_waves = sun[1].data['WAVELENGTH'] * u.Angstrom 
        # sfd = spectral flux density
        sun_sfd = sun[1].data['FLUX'] * u.erg * u.second**-1 * u.cm**-2 * u.Angstrom**-1

        self.waves = sun_waves.to(u.micron)
        
        # Covert to zodiacal light spectrym by following the normalisation and reddening
        # prescription of Leinert et al (1997) with the revised parameters from 
        # Aldering (2001), as used in the HST ETC (Giavalsico, Sahi, Bohlin (2202)).
        
        # Reddening factor
        rfactor = np.where(sun_waves < ZodiacalLight.lambda_c, \
                           1.0 + ZodiacalLight.f_blue * np.log(sun_waves/ZodiacalLight.lambda_c), \
                           1.0 + ZodiacalLight.f_red * np.log(sun_waves/ZodiacalLight.lambda_c))
        # Apply normalisation and reddening
        sfd = sun_sfd * ZodiacalLight.zl_normalisation * rfactor
        # #DownWithErgs
        self.sfd = sfd.to(u.Watt * u.m**-2 * u.arcsecond**-2 * u.micron**-1)
        # Also calculate in photon spectral flux density units. Fudge needed because equivalencies
        # don't currently include surface brightness units.
        fudge = sfd * u.arcsecond**2
        fudge = fudge.to(u.photon * u.second**-1 * u.m**-2 *  u.micron**-1, equivalencies=u.spectral_density(self.waves))
        self.photon_sfd = fudge / u.arcsecond**2
        
