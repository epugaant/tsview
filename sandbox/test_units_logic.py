import astropy.units as u
import numpy as np
from synphot import units, SourceSpectrum

flux = 203655.62323629577 #electrons/s
flux_error = 76.00599243480103 # electrons/s
mag_vegamag = 12.415125851786478

zeropt = {'G': {'VEGAMAG': {'zp': 25.6873668671 * units.VEGAMAG, 'e_zP': 0.0027553202* units.VEGAMAG, 'lamb': 6251.50 * u.AA, 'f_zp_nu': 3228.75 * u.Jy}, 'AB': {'zp': 25.8010446445 * u.ABmag, 'e_zP': 0.0027590522 * u.ABmag, 'lamb': 6251.50 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
          'BP': {'VEGAMAG': {'zp': 25.3385422158 * units.VEGAMAG, 'e_zP': 0.0027901700* units.VEGAMAG, 'lamb': 5124.20 * u.AA, 'f_zp_nu': 3552.01 * u.Jy}, 'AB': {'zp': 25.3539555559 * u.ABmag, 'e_zP': 0.0023065687 * u.ABmag, 'lamb': 5124.20 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
          'RP': {'VEGAMAG': {'zp': 24.7478955012 * units.VEGAMAG, 'e_zP': 0.0037793818* units.VEGAMAG, 'lamb': 7829.65 * u.AA, 'f_zp_nu': 2554.95 * u.Jy}, 'AB': {'zp': 25.1039837393 * u.ABmag, 'e_zP': 0.0015800349 * u.ABmag, 'lamb': 7829.65 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
          }
_u_vega = u.def_unit('VEGA')
VEGAMAG = u.mag(_u_vega)
ABMAG = u.ABmag
# Register with astropy units
u.add_enabled_units([VEGAMAG, ABMAG])

G_AB = u.Magnitude(flux* u.electron /u.s) + zeropt['G']['AB']['zp']*u.ABmag
G_err_AB = u.Magnitude(flux_error* u.electron /u.s) + zeropt['G']['AB']['zp']*u.ABmag

F_G_mJy = (G_AB.value * u.ABmag).to(u.mJy)
F_err_G_mJy = (G_err_AB.value * u.ABmag).to(u.mJy)

original_unit = u.Unit('electron/s')
zp_unit = u.ABmag
target_unit = u.mJy

#check zp_unit 
VEGA = u.def_unit('mag ({0})'.format('VEGA'), 1 * units.VEGAMAG)
u.add_enabled_units([VEGA])
isinstance((1* u.Unit('mag ({0})'.format('AB'))), u.Magnitude) # True
isinstance((1* units.VEGAMAG), u.Magnitude) # True


eqv = u.spectral_density(zeropt['G']['AB']['lamb_AA']* u.AA)

(values * u.Unit(original_units)).to_value(u.Unit(target_units), equivalencies=eqv)


_wave = [4956.8, 4959.55, 4962.3] * u.AA
_flux_photlam = [9.7654e-3, 1.003896e-2, 9.78473e-3] * units.PHOTLAM
_flux_photnu = [8.00335589e-14, 8.23668949e-14, 8.03700310e-14] * units.PHOTNU
_flux_flam = [3.9135e-14, 4.0209e-14, 3.9169e-14] * units.FLAM
_flux_fnu = [3.20735792e-25, 3.29903646e-25, 3.21727226e-25] * units.FNU
_flux_jy = [3.20735792e-2, 3.29903646e-2, 3.21727226e-2] * u.Jy
_flux_count = [1214.88479883, 1248.91795446, 1217.28946691] * u.count
_flux_stmag = [12.41858665, 12.38919182, 12.41764379] * u.STmag
_flux_abmag = [12.63463143, 12.60403221, 12.63128047] * u.ABmag
_flux_obmag = [-7.71133775, -7.74133477, -7.71348466] * units.OBMAG
_flux_vegamag = [12.72810665, 12.69861694, 12.72605148] * units.VEGAMAG

#('in_q', 'out_u', 'ans', 'support_scalar')
(_flux_abmag, units.PHOTLAM, _flux_photlam, True)
result = units.convert_flux(_wave, _flux_abmag, u.Jy)
vega = SourceSpectrum.from_vega()  # For unit conversion  
result = units.convert_flux(_wave, _flux_vegamag, u.Jy, vegaspec=vega)

import astropy.units as u
import numpy as np
from synphot import units, SourceSpectrum, SpectralElement

units.convert_flux([5124.2 , 6251.5 , 7829.65] * u.AA, [0, 0, 0] * units.VEGAMAG, u.Jy, vegaspec=vega)

import matplotlib.pyplot as plt
from synphot import SourceSpectrum, SpectralElement, Observation, units
from synphot.specio import read_remote_spec
from astropy.io import ascii
from specutils import Spectrum1D
import os


path = '/Users/epuga/ESDC/Gaia/Photometry/GaiaEDR3_passbands_zeropoints_version2'
table = ascii.read(os.path.join(path,"passband.dat"), readme=os.path.join(path, 'ReadMe'))
table.meta['expr'] = filename = 'gaia_G.fits'
table['WAVELENGTH'] = table['lambda']
table['THROUGHPUT'] = table['GPb']
table['THROUGHPUT'].unit = ''
table['WAVELENGTH', 'THROUGHPUT'].write(os.path.join(path,filename), format='fits', overwrite=True) 
table.remove_columns(['THROUGHPUT'])
table.meta['expr'] = filename = 'gaia_BP.fits'
table['THROUGHPUT'] = table['BPPb']
table['WAVELENGTH', 'THROUGHPUT'].write(os.path.join(path,filename), format='fits') 
table.meta['expr'] = filename = 'gaia_RP.fits'
table['THROUGHPUT'] = table['RPPb']
table['WAVELENGTH', 'THROUGHPUT'].write(os.path.join(path,filename), format='fits') 


vega_internal = True
norm = True
bp_v = SpectralElement.from_filter('johnson_v')  
if vega_internal:
    vega = SourceSpectrum.from_vega()  # For unit conversion  
else:
    hdr, x, f = read_remote_spec('https://ssb.stsci.edu/cdbs/calspec/alpha_lyr_mod_004.fits', encoding='binary')
    spec = Spectrum1D(spectral_axis=x, flux=f)
    vega = SourceSpectrum.from_spectrum1d(spec)
bp_g = SpectralElement.from_file('/Users/epuga/ESDC/Gaia/Photometry/GAIA_GAIA3.Gbp.dat')
sp_norm = vega.normalize(0.0 * units.VEGAMAG, bp_v, vegaspec=vega)  
if norm:
    obs = Observation(sp_norm, bp_g) 
else:
    obs = Observation(vega, bp_g) 
#This is the function according to https://core.ac.uk/download/pdf/156679601.pdf
print(bp_g.pivot(), obs.effstim(flux_unit='jy'))
obs.effstim(flux_unit='flam')

val = obs.integrate(flux_unit='FLAM', binset=range(3200,10500))/bp_g.integrate()
print(bp_g.pivot(), units.convert_flux(bp_g.pivot(), val, 'jy'))

#backward check for normalization
obs_v = Observation(sp_norm, bp_v)
obs_v.effstim(flux_unit='mag(VEGA)', vegaspec=sp_norm)

# G band
sp_norm(6251.5 * u.AA, flux_unit=u.Jy) 
wave = sp_norm.waveset
plt.plot(wave, vega(wave, flux_unit='jy'), 'k', wave, sp_norm(wave, flux_unit='jy'), 'b')  
plt.plot(wave, vega(wave, flux_unit='jy'), 'r')
plt.xlim(1000, 30000)  
plt.xlabel('Wavelength (Angstrom)')  
plt.ylabel('Flux (Jy)')  
#plt.title(vega.meta['expr'])  
plt.legend(['Original', 'Normalized'], loc='upper right')  
plt.show()

#BP, G, RP
units.convert_flux([5124.2 , 6251.5 , 7829.65] * u.AA, [0, 0, 0] * units.VEGAMAG, u.Jy, vegaspec=sp_norm)

plt.plot(wave, sp_norm(wave, flux_unit='jy'), 'k')
plt.xlim(1000, 30000)
plt.xlabel('Wavelength (Angstrom)')
plt.ylabel('Flux (Jy)')
val = obs.integrate(flux_unit='FLAM', binset=range(3200,10500))/bp_g.integrate()
plt.plot(bp_g.pivot(), units.convert_flux(bp_g.pivot(), val, 'jy'), 'co')
fnu_0 = units.convert_flux([5109.71, 6217.51, 7769.02] * u.AA, [0, 0, 0] * units.VEGAMAG, u.Jy, vegaspec=sp_norm)
plt.plot([5109.71, 6217.51, 7769.02] * u.AA, fnu_0, 'b*')
'''This is the normalization that convert_flux uses, but it is not the zeropoint calculated'''

#f_zp_nu = column_from_dict_index(t, 'band', gaia, '**.{{}}.**.{0}.**.{1}'.format('VEGAMAG', 'f_zp_nu'))
f_zp_nu = [3552.01, 3228.75, 2554.95] * u.Jy
plt.plot([5109.71, 6217.51, 7769.02] * u.AA, f_zp_nu, 'y+')
plt.show()

#plot of 
plt.plot(sp_norm.waveset, sp_norm(sp_norm.waveset, flux_unit='jy'), 'k')
norm = units.convert_flux([6217.51] * u.AA, [0] * units.VEGAMAG, u.Jy, vegaspec=sp_norm)
plt.plot(bp_g.waveset, norm*bp_g(bp_g.waveset), 'r')
plt.plot(obs.waveset, obs(obs.waveset, flux_unit='jy'), 'b')
plt.xlim(3000, 15000)
plt.show()