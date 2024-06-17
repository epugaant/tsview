#------Physical units-------

u.ABflux.physical_type
#PhysicalType('spectral flux density')

u.ABflux.bases
#[Unit("AB")]

u.ABflux.represents
#Unit("3.63078e-20 erg / (Hz s cm2)")

u.ABflux.decompose()
#Unit("3.63078e-23 kg / s2")

u.ABflux.si
#Unit("3.63078e-23 kg / s2")

u.ABflux.cgs
#Unit("3.63078e-20 g / s2")

u.ABflux.in_units(u.Jy)
#3630.7805477010033

u.ABflux.in_units(u.Watt/(u.m*u.m*u.Hz))
#3.6307805477010035e-23

u.ABflux.to(u.Jy)
#3630.7805477010033


#-------Logarithmic Unit---------

u.ABmag.physical_unit
#Unit("AB")

u.ABmag.physical_type
#PhysicalType('spectral flux density')

u.ABmag.decompose()
#Unit("mag(3.63078e-23 kg / s2)")

u.ABflux.is_equivalent(u.Jy)
#True

u.ABflux.is_equivalent(u.mag)
#False

u.ABmag.equivalencies
#[(Unit("mag(AB)"),
#  Unit("AB"),
#  <bound method LogUnit.to_physical of Unit("mag(AB)")>,
#  <bound method LogUnit.from_physical of Unit("mag(AB)")>)]

u.ABmag.from_physical(2000 )
#-8.252574989159953

#Wich is equivalent to doing
u.Magnitude(2000 *u.ct/u.s)
#<Magnitude -8.25257499 mag(ct / s)>

u.ABmag.to_physical(-8.25257499 )
#2000.0000015474238






# Tests to learn astropy units Photometry (Magnitudes and Other Logarithmic Units)
        
'''The zero_point_flux() equivalency provides a way to move between photometric systems 
(i.e., those defined relative to a particular zero-point flux) and absolute fluxes.'''
from astropy.units import u

u.Jy.physical_type
u.Jy.find_equivalent_units()

'''Maggies - a linear flux unit that is the flux for a mag=0 object. 
To tie this onto a specific calibrated unit system, the zero_point_flux equivalency should be used.'''

'''An equivalency for converting linear flux units ("maggys") defined 
relative to a standard source into a standardized system.'''
u.zero_point_flux(3631.1 * u.Jy)
#[(Unit("mgy"), Unit("3631.1 Jy"))]


# from photometry test: 
(1 * u.nmgy).to(u.ABflux, u.zero_point_flux(1 * u.ABflux))
#<Quantity 1.e-09 AB> so this is f0/f but is not the 3631e-9 from the example
u.Magnitude((1 * u.nmgy).to(u.ABflux, u.zero_point_flux(1 * u.ABflux))) 
#<Magnitude 22.5 mag(AB)>

# From equivalencies, this is the use case they give 
u.Magnitude((1.2 * u.nanomaggy).to(u.ABflux, u.zero_point_flux(3631.1 * u.Jy))) 
#<Magnitude 22.30195136 mag(AB)>

 # Definition of zeropoint flux
 #flux zero-point in AB system is 3631 Jy; 1Jy=10E-26 W/Hz/m**2
 
'''Another way to express these zeropoints is to say that an object 
with f = 3.63 10-20 erg cm-2 s-1 Hz-1 will have mAB=0 in every filter, 
and an object with f = 3.63 10-9 erg cm-2 s-1 A-1 will have mST=0 in every filter.'''
 
# STmag (lambda)
(0. * u.STmag).to(u.erg/u.s/u.cm**2/u.AA)  
#<Quantity 3.63078055e-09 erg / (Angstrom s cm2)>
(-21.1 * u.STmag).to(u.erg/u.s/u.cm**2/u.AA) 
#<Quantity 1. erg / (Angstrom s cm2)>

# ABmag (nu)
(-48.60 * u.ABmag).to(u.erg/u.s/u.cm**2/u.Hz) 
#<Quantity 1. erg / (Hz s cm2)>
(0. * u.ABmag).to(u.erg/u.s/u.cm**2/u.Hz) 
#<Quantity 3.63078055e-20 erg / (Hz s cm2)>

#ABflux in Jy to mag
u.Magnitude(3631*u.Jy)
#<Magnitude -8.90006562 mag(Jy)>

#ABflux in W/Hz/m^2 to mag
u.Magnitude(3631e-26*u.Watt/u.Hz/u.m**2)
#<Magnitude 56.09993438 mag(W / (Hz m2))>

#maggys in Jy to mag. It differs from the one in the University of wyoming presentation
u.Magnitude(1e9*u.Jy)
#<Magnitude -22.5 mag(Jy)>

'''The relationship between ABMAG and STMAG is:

    ABMAG = STMAG - 5 log (PHOTPLAM) + 18.692'''

'''PHOTFLAM is the flux of a source with constant flux per unit wavelength (in erg s-1 cm-2 �-1) which produces a count rate of 1 DN per second.'''
'''PHOTPLAM is the pivot wavelength in Angstroms, where pivot wavelength is a measure of the effective wavelength of a filter'''

# u.ABmag (nu) --> F(lambda) [u.erg/u.s/u.cm**2/u.AA] 
#astropy slack example
flam = (20 * u.ABmag).to(u.erg/u.s/u.cm**2/u.AA, u.spectral_density(8000 * u.AA))

import astropy.units as u
from astroquery.simbad import Simbad
from astropy.table import Table

filename= '/Users/epuga/ESDC/TSViz/data/astropy_units/lwr05639mxlo_vo.fits'
t_lwr = Table.read(filename)
print(t_lwr)

wav_UV = t_lwr['WAVE'][0,].quantity #comes packed in two
UVflux = t_lwr['FLUX'][0,].quantity

custom_query = Simbad()
custom_query.add_votable_fields('fluxdata(U)','fluxdata(B)','fluxdata(V)')
phot_table=custom_query.query_object('HD 147933')
Umag=phot_table['FLUX_U']
Bmag=phot_table['FLUX_B']
Vmag=phot_table['FLUX_V']

#Vega system is default from http://ned.ipac.caltech.edu/help/photoband.lst
wav_U = 0.3660 * u.micron 
zeroflux_U_nu = 1.81E-23 * u.Watt/(u.m*u.m*u.Hz)
wav_B = 0.4400 * u.micron
zeroflux_B_nu = 4.26E-23 * u.Watt/(u.m*u.m*u.Hz)
wav_V = 0.5530 * u.micron
zeroflux_V_nu = 3.64E-23 * u.Watt/(u.m*u.m*u.Hz)

zeroflux_U = zeroflux_U_nu.to(u.erg/u.AA/u.cm/u.cm/u.s, 
                              equivalencies=u.spectral_density(wav_U))
zeroflux_B = zeroflux_B_nu.to(u.erg/u.AA/u.cm/u.cm/u.s, 
                              equivalencies=u.spectral_density(wav_B))
zeroflux_V = zeroflux_V_nu.to(u.erg/u.AA/u.cm/u.cm/u.s, 
                              equivalencies=u.spectral_density(wav_V))

Uflux = zeroflux_U * 10.**(-0.4*Umag)
Bflux = zeroflux_B * 10.**(-0.4*Bmag)
Vflux = zeroflux_V * 10.**(-0.4*Vmag)

#For Vega Magnitudes
Umag = 4.3 
flam_zero_U = (1.81E-23 * u.Watt/(u.m*u.m*u.Hz)).to(u.erg/u.AA/u.cm/u.cm/u.s, 
                              equivalencies=u.spectral_density(wav_U)) 
flam_U = flam_zero_U * 10.**(-0.4*Umag)
#<Quantity [7.718572e-11] erg / (Angstrom s cm2)>

(Umag.value * u.mag).to(u.erg/u.s/u.cm**2/u.AA, u.zero_point_flux((1.81E-23 * u.Watt/(u.m*u.m*u.Hz)).to(u.erg/u.AA/u.cm/u.cm/u.s, 
                              equivalencies=u.spectral_density(wav_U))))

(Umag.value * u.mag).to(u.erg/u.AA/u.cm/u.cm/u.s, u.spectral_density(wav_U))

from synphot import units
Umag.value *units.VEGAMAG
(Umag.value * units.VEGAMAG).to(u.erg/u.s/u.cm**2/u.AA, u.spectral_density(wav_U)) #ERROR

#Following the example of test_photometric
ST_base_unit = u.erg * u.cm**-2 / u.s / u.AA
(10 * u.mgy).to(u.STflux, u.zero_point_flux(1 * ST_base_unit))
#this is not the 10 * ST_base_unit

AB_base_unit = u.erg/u.s/u.cm**2/u.AA
(Umag.value * u.ABmag).to(u.ABflux, u.zero_point_flux(1 * AB_base_unit))