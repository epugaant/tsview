import astropy.units as u

tint = 1000.*u.s
# count rates for three sources in filters b and v
cr_b = ([3000., 100., 15.] * u.ct) / tint
cr_v = ([4000., 90., 25.] * u.ct) / tint
b_i, v_i = u.Magnitude(cr_b), u.Magnitude(cr_v)
#instrument magnitudes
b_i, v_i  

# first star has known magnitudes (standard star)
b_ref, v_ref = 17.2 * u.STmag, 17.0 * u.STmag
b_ref, v_ref  

#consider the fist star the zeropoint 
zp_b, zp_v = b_ref - b_i[0], v_ref - v_i[0]
zp_b, zp_v  
# (<Magnitude 18.392803...T s / ct)>, <Magnitude 18.505149...T s / ct)>)

#absolute magnitudes in ST standard system
B, V = b_i + zp_b, v_i + zp_v
B, V  
V.to(u.ABmag, u.spectral_density(5500.*u.AA))  

#STmag to ABmag
V.to(u.ABmag, u.spectral_density(5500.*u.AA))

flam = V.to(u.erg/u.s/u.cm**2/u.AA)
flam 

lam = 5500 * u.AA
fnu = V.to(u.erg/u.s/u.cm**2/u.Hz, u.spectral_density(lam))
fnu



import numpy as np

flux = 203655.62323629577 #electrons/s
flux_error = 76.00599243480103 # electrons/s
mag_vegamag = 12.415125851786478


# generate magerr from fluxerr and flux
magerr = 2.5 / np.log(10) * flux_error / flux

# compute flux and flux error in mJy for what they are using, thinking that mag_vegamag is in AB system
flux_mJy = 10 ** (-0.4 * (mag_vegamag - 23.9)) / 1e3  # in mJy # 39.25994250790246
fluxerr_mJy = magerr / 2.5 * np.log(10) * flux_mJy  # in mJy


zeropt = {'G': {'VEGAMAG': {'zp': 25.6873668671, 'e_zP': 0.0027553202, 'lamb_AA': 6251.50, 'f_zp_nu_Jy': 3599.12}, 'AB': {'zp': 25.8010446445, 'e_zP': 0.0027590522, 'lamb_AA': 6251.50, 'f_zp_nu_Jy': 3631}}, 
          'BP': {'VEGAMAG': {'zp': 25.3385422158, 'e_zP': 0.0027901700, 'lamb_AA': 5124.20, 'f_zp_nu_Jy': 3256.01}, 'AB': {'zp': 25.3539555559, 'e_zP': 0.0023065687, 'lamb_AA': 5124.20, 'f_zp_nu_Jy': 3631}}, 
          'RP': {'VEGAMAG': {'zp': 24.7478955012, 'e_zP': 0.0037793818, 'lamb_AA': 7829.65, 'f_zp_nu_Jy': 2530.35}, 'AB': {'zp': 25.1039837393, 'e_zP': 0.0015800349, 'lamb_AA': 7829.65, 'f_zp_nu_Jy': 3631}}, 
          }
G_vegamag = u.Magnitude(flux* u.electron /u.s) + zeropt['G']['VEGAMAG']['zp']*u.mag
G_err_vegamag = u.Magnitude(flux_error* u.electron /u.s) + zeropt['G']['VEGAMAG']['zp']*u.mag
G_AB = u.Magnitude(flux* u.electron /u.s) + zeropt['G']['AB']['zp']*u.ABmag

F_G_mJy_hector = zeropt['G']['VEGAMAG']['f_zp_nu_Jy']*10.**(-0.4*G_vegamag.value)*1E3 # test to check that 38.91759427948141 it is equivalent to flux_mJy 

flux_mJy = 10 ** (-0.4 * (G_AB.value - 23.9)) / 1e3
f = u.Dex(-0.4*(mag_vegamag - 23.9))
f.physical/1e3
# let's try writing it in astropy.units
F_G_mJy = (G_AB.value * u.ABmag).to(u.mJy)

F_G_mJy = (G_vegamag.value * u.ABmag).to(u.mJy)
F_err_G_mJy = (G_err_vegamag.value * u.ABmag).to(u.mJy)


