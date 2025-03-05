import numpy as np
import scipy.constants
from dataclasses import dataclass

#------------------------------------------------------------------------
#------------------------------------------------------------------------
@dataclass
class spec:
    a: float # aparture area [cm2]
    pix_sz: float # detector pixel size [micro-m]
    d: float    # dispersion [A/pix]
    m: float    # plate scale [arcsec/pix]
    w: float    # slit width [arcsec]
    n_w: int   # number of pixel in dispersion direction
    n_s: int   # number of pixel in spatial direction
    cr_d_min: float # dark current count rate [/s/pix]
    cr_d_max: float
    cr_r_min: float # penetrating high energy particle count rate [/s/pix]
    cr_r_max: float
    description: str
    I_Lya: float  #  Typical geocorona brightness * Ly-a
    I_1304: float  #  Typical geocorona brightness * OI 1304
    I_1356: float  #  Typical geocorona brightness * OI 1356

#------------------------------------------------------------------------
#------------------------------------------------------------------------
# LAPYUTA spec
def get_lap_spec(inst='mrs', model='baseline'):

    # Typical geocorona brightness for LAPYUTA [erg/s/cm2/arcsec2/A]
    I_Lya = 6.1e-14   # 2kR (HST daytime / 10)
    I_1304 = 5.7e-17  # 2R  (HISAKI)
    I_1356 = 0.0

    a = (60.0/2.0)**2 * scipy.constants.pi

    if inst == 'mrs':
        pix_sz = 15.0  # detector pixel size [micro-meter/pixel]
        d = 0.1 
        m = 0.055 
        w = 1.0 
        n_w = 10                 # 1A幅 (1/d)
        n_s = 18                 # 1arcsec幅 (1/m)
        description = 'LAPYUTA MRS ' + model
    elif inst == 'hrs':
        pix_sz = 15.0  # detector pixel size [micro-meter/pixel]
        d = 0.015
        m = 0.3 
        w = 1.0 
        n_w = 67                 # 1A幅 (1/d)
        n_s = 3                  # 1arcsec幅 (1/m)
        description = 'LAPYUTA HRS ' + model

    # detector noise
    cr_r_min = 0.3  * (pix_sz * 1e-4)**2  # /s/cm2 -> /s/pixel
    cr_r_max = 10.0 * (pix_sz * 1e-4)**2
    cr_d_min = 0.43 * (pix_sz * 1e-4)**2  # /s/cm2 -> /s/pixel
    cr_d_max = 0.43 * (pix_sz * 1e-4)**2

    lap = spec(
        a = a,
        pix_sz = pix_sz,
        d = d, 
        m = m, 
        w = w, 
        n_w = n_w,
        n_s = n_s, 
        cr_r_min = cr_r_min,
        cr_r_max = cr_r_max,
        cr_d_min = cr_d_min,
        cr_d_max = cr_d_max,
        description = description,
        I_Lya = I_Lya,
        I_1304 = I_1304,
        I_1356 = I_1356
    )

    return lap

#------------------------------------------------------------------------
# HST spec
def get_hst_spec(inst='STIS_G140M', cycle=24):

    # Typical geocorona brightness for HST [erg/s/cm2/arcsec2/A]
    I_Lya = 6.1e-13   # 20kR (HST daytime)
    I_1304 = 5.7e-14 # 2kR  (HST daytime)
    I_1356 = 5.0e-15 # 0.2kR(HST daytime)

    if inst == 'STIS_G140L':
        # STIS G140L
        a = (240.0/2.0)**2 * scipy.constants.pi
        pix_sz = 25.0 # FUV-MAMA
        d = 0.6
        m = 0.0246
        w = 1.0
        n_w = 2                # 1A幅
        n_s = 2
        description = 'HST STIS G140L'
    else:
        # STIS G140M
        a = (240.0/2.0)**2 * scipy.constants.pi
        pix_sz = 25.0 # FUV-MAMA
        d = 0.05
        m = 0.029
        w = 1.0
        n_w = 20               # 1A幅 (1/d)
        n_s = 34               # 1arcsec幅 (1/m)
        description = 'HST STIS G140M'

    # detector noise level
    if cycle == 8:
        cr_d_min = 5e-6        # STIS IHB C8 (1998)
        cr_d_max = 1e-5        # STIS IHB C8 (1998)
        cr_r_min = 0.0
        cr_r_max = 0.0
    else:
        cr_d_min = 7e-6        # STIS IHB C24 (2024)
        cr_d_max = 6e-4        # STIS IHB C24 (2024)
        cr_r_min = 0.0
        cr_r_max = 0.0

    hst = spec(
        a = a,
        pix_sz = pix_sz,
        d = d, 
        m = m, 
        w = w, 
        n_w = n_w,
        n_s = n_s, 
        cr_r_min = cr_r_min,
        cr_r_max = cr_r_max,
        cr_d_min = cr_d_min,
        cr_d_max = cr_d_max,
        description = description,
        I_Lya = I_Lya,
        I_1304 = I_1304,
        I_1356 = I_1356
        )

    return hst

#------------------------------------------------------------------------
# UVEX spec
#
# Detector noise
# Delta-doped Silicon UV/Visible Detectors at JPL
# ex) Delta-doped electron-multiplying CCDs for FIREBall-2 (Gillian Kyne et al., J. Astron. Telesc. Instrum. Syst. 6(1), 011007 (2020), doi: 10.1117/1.JATIS.6.1.011007)
#  Pixel size: 13um
#　CIC  = 10^-3 [electron/pixel/frame]  (clock induced count)
#　Dark = 10^-2 [electron/pixel/hour]
#
# The Optical Design for the Ultraviolet Explorer (UVEX) Mission: a next-generation wide-field UV telecope by Jason Fucik (SPIE2024, https://doi.org/10.1117/12.3020427)
#  Pixel size: 10um, Slit Width: 2", 4", 8", Linear Dispersion 0.4 A/pix, Plate Scale 0.9"/pix (Table 1)
#
def get_uvex_spec():

    # Typical geocorona brightness for UVEX [erg/s/cm2/arcsec2/A]
    I_Lya = 6.1e-14/4.0   # 0.5kR (HST daytime / 40)
    I_1304 = 0.0
    I_1356 = 0.0

    a = (75.0/2.0)**2 * scipy.constants.pi
    pix_sz = 13.0  # detector pixel size [micro-meter/pixel]
    d = 0.4              # A/pix
    m = 0.9              # arcsec/pix
    w = 1.0              # arcsec
#    w = [2.0, 4.0, 8.0]  # arcsec
    n_w = 2              # 1A幅 (1/d)
    n_s = 1              # 1arcsec幅 (1/m)
    cr_r_min = 0.0
    cr_r_max = 0.0
    cr_d_min = 1e-3*(10.0/13.0)**2 # /s/pixel
    cr_d_max = 1e-3*(10.0/13.0)**2
    description = 'UVEX LSS'

    uvex = spec(
        a = a,
        pix_sz = pix_sz,
        d = d, 
        m = m, 
        w = w, 
        n_w = n_w,
        n_s = n_s, 
        cr_r_min = cr_r_min,
        cr_r_max = cr_r_max,
        cr_d_min = cr_d_min,
        cr_d_max = cr_d_max,
        description = description,
        I_Lya = I_Lya,
        I_1304 = I_1304,
        I_1356 = I_1356
        )

    return uvex

#------------------------------------------------------------------------
#------------------------------------------------------------------------
# LAPYUTA Ae
def get_lap_ae(inst='mrs', model='baseline'):
    filename = 'LAPYUTA/LAPYUTA_sensitivity_20240328.csv'
    ae_array = np.loadtxt(filename, skiprows=1, delimiter=',')
    wl_ = ae_array[:,0] * 10.0  # [A]

    d_wl = 0.005  # [A]
    nw = int((wl_[-1] - wl_[0]) / d_wl)
    lap_wl = np.linspace(wl_[0], wl_[-1], nw)
    if model == 'baseline':
        lap_ae = scipy.interpolate.interp1d(wl_,ae_array[:,1])(lap_wl) # MRS Baseline
    else: # target
        lap_ae = scipy.interpolate.interp1d(wl_, ae_array[:,2])(lap_wl) # MRS Target
    #ae_mrs_baseline = scipy.interpolate.interp1d(wl_,ae_array[:,1])(lap_wl)
    #ae_mrs_target   = scipy.interpolate.interp1d(wl_,ae_array[:,2])(lap_wl)
    #ae_uvi_baseline = scipy.interpolate.interp1d(wl_,ae_array[:,3])(lap_wl)
    #ae_uvi_target   = scipy.interpolate.interp1d(wl_,ae_array[:,4])(lap_wl)

    if inst == 'hrs':
        lap_ae = lap_ae * 0.5

    return lap_wl, lap_ae

#------------------------------------------------------------------------
# HST Ae
def get_hst_ae(inst='STIS_G140M'):

    filename = 'HST/STIS-G140M_Throughput.csv'
    ae_array = np.loadtxt(filename, skiprows=0, delimiter=',')
    wl_ = ae_array[:,0]  # [A]
    tp  = ae_array[:,1] * 0.01 # Throughput

    d_wl = 0.005  # [A]
    nw = int((wl_[-1] - wl_[0]) / d_wl)
    hst_wl = np.linspace(wl_[0], wl_[-1], nw)

    a = (240.0/2.0)**2 * scipy.constants.pi
    hst_ae = scipy.interpolate.interp1d(wl_, tp)(hst_wl) * a

    return hst_wl, hst_ae

#------------------------------------------------------------------------
# UVEX Ae
# https://uvex2023.caltech.edu/system/media_files/binaries/6/original/UVEX_slides_bg_-_Brian_Grefenstette.pdf?1680041999
def get_uvex_ae():

    filename = 'UVEX/UVEX_LSS_Brian_2023.csv'
    ae_array = np.loadtxt(filename, skiprows=2, delimiter=',')
    wl_ = ae_array[:,0]  # [A]
    tp  = ae_array[:,1]  # Throughput

    d_wl = 0.005  # [A]
    nw = int((wl_[-1] - wl_[0]) / d_wl)
    uvex_wl = np.linspace(wl_[0], wl_[-1], nw)

    a = (75.0/2.0)**2 * scipy.constants.pi
    uvex_ae = scipy.interpolate.interp1d(wl_, tp)(uvex_wl) * a

    return uvex_wl, uvex_ae

#------------------------------------------------------------------------
# HST sensitivity
def get_hst_sensitivity(inst='STIS_G140M'):

    # HST (STIS Instrument Handbook for Cycle 32)
    filename = 'HST/HST_STIS_G140M.csv'
    hst_array = np.loadtxt(filename, skiprows=1, delimiter=',')
    wl = hst_array[:,0]

    d_wl = 0.005  # [A]
    nw = int((wl[-1] - wl[0]) / d_wl)
    hst_wl = np.linspace(wl[0], wl[-1], nw)
    hst_sp = scipy.interpolate.interp1d(wl, hst_array[:,-1])(hst_wl)

    return hst_wl, hst_sp

#------------------------------------------------------------------------
#------------------------------------------------------------------------
def line_profile(amplitude, center, fwhm, x):
    sgm = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    arg = -(x-center)**2 / (2.0 * sgm**2)
    y = amplitude * np.exp(arg) / sgm / np.sqrt(2.0 * scipy.constants.pi)
    return y

#------------------------------------------------------------------------
def get_geocorona(I_Lya, I_1304, I_1356, wl):

    I  = line_profile(I_Lya, 1215.67, 0.04, wl) 
    I += line_profile(I_1304*3.0/6.0, 1302.2, 0.013, wl)
    I += line_profile(I_1304*2.0/6.0, 1304.9, 0.013, wl)
    I += line_profile(I_1304*1.0/6.0, 1306.0, 0.013, wl)
    I += line_profile(I_1356,         1356.0, 0.013, wl)

    return I

#------------------------------------------------------------------------
def get_sensitivity(spec, ae, wl, resolved=1):
    h = scipy.constants.h * 1e7  # J s -> erg s
    c = scipy.constants.c * 1e2  # m/s -> cm/s

    if resolved == 1:   # resolved line or continnum
        # point source sensitivity [(counts/s/pix) / (erg/cm2/s/A)]
        sp = ae * spec.d * (wl * 1e-8)/(h * c)
        #   cm^2   A/pix     A->cm      erg cm

        # diffuse source sensitivity [(counts/s/pix2) / (erg/cm2/s/A/arcsec2)]
        sd = sp * spec.m * spec.w
        #      [arcsec/pix] [arcsec]
    else:                # un-resolved line
        # point source sensitivity [(counts/s/pix) / (erg/cm2/s)]
        sp = ae * spec.d * (wl * 1e-8)/(h * c)
        #   cm^2   A/pix     A->cm      erg cm

        # diffuse source sensitivity [(counts/s/pix2) / (erg/cm2/s/arcsec2)]
        sd = ae * (wl * 1e-8)/(h * c) * spec.m * spec.m
        #      [arcsec/pix] [arcsec]

    return sp, sd
