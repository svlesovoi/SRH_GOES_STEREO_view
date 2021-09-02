#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 27 17:10:14 2021

@author: sergeyvlesovoi
"""


import netCDF4 as nc
import numpy as NP
import cftime
import pylab as PL
from datetime import datetime
from astropy.io import fits
from optparse import OptionParser
import os
import urllib.request
from scipy.io import readsav

def hm_format(t, pos):
  hh = int(t / 60.)
  t -= hh*60.
  mm = int(t)
  return '%02d:%02d' % (hh,mm);

def freq_format(f, pos):
    if (f < 368):
        return '%3.1f' % (stereoA['frequencies'][int(f)] * 1e-3)
    else:
        return ''

def srhTime2HHMMSS(srhTimeSecs):
    HH = int(srhTimeSecs // 3600)
    srhTimeSecs -= HH*3600
    MM = int(srhTimeSecs // 60)
    srhTimeSecs -= MM*60
    SS = int(srhTimeSecs)
    return [HH,MM,SS]

def smooth(x, length):
    kern = NP.ones(length)/length
    x_smooth = NP.convolve(x, kern, mode='same')
    return x_smooth

def GOESnameFromDate(requestDate):
    GOESprefix = 'sci_xrsf-l2-flx1s_g16_d'
    GOESsuffix = '_v2-1-0.nc'
    return GOESprefix + requestDate.replace('/','') + GOESsuffix

def SRHnameFromDate(requestDate):
    SRHprefix = 'srh_cp_'
    SRHsuffix = '.fits'
    return SRHprefix + requestDate.replace('/','') + SRHsuffix
    
smoothWindow = 21

parser = OptionParser()
parser.add_option("-d", "--date", dest="requestDate", default = '2021/08/05')
(goes_srh_options, goes_srh_args) = parser.parse_args()

requestDate = goes_srh_options.requestDate
srhCpFitPath = SRHnameFromDate(requestDate) #goes_srh_options.cpFitPath
goesPath =  GOESnameFromDate(requestDate) #goes_srh_options.goesPath

if not os.path.exists(srhCpFitPath):
    srhURL = 'http://archive.rao.istp.ac.ru/SRH/' + srhCpFitPath
    urllib.request.urlretrieve(srhURL, srhCpFitPath)

stereoPath = 'swaves_average_' + requestDate.replace('/','') + '_a.sav'
if not os.path.exists(stereoPath):
    stereoURL = 'https://solar-radio.gsfc.nasa.gov/data/stereo/new_summary/' + requestDate[0:4] + '/' + stereoPath
    urllib.request.urlretrieve(stereoURL, stereoPath)

stereoA = readsav(stereoPath)

if not os.path.exists(goesPath):
    dateList = requestDate.split('/')
    url_path = 'https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/goes16/l2/data/xrsf-l2-flx1s_science/'
    urllib.request.urlretrieve(url_path + dateList[0] + '/' + dateList[1] +'/' + goesPath, goesPath)

goesFile = nc.Dataset(goesPath)
goesDatetime = cftime.num2pydate(goesFile.variables["time"][:], goesFile["time"].units)
platform = getattr(goesFile, "platform")
var_name = ["xrsa_flux", "xrsb_flux"]

srhCpFitFile = fits.open(srhCpFitPath)

dateList = srhCpFitFile[0].header['DATE-OBS'].split('-')
srhFreqList = srhCpFitFile[1].data['frequencies']
srhTime = srhCpFitFile[2].data['time']
srhCorrI = srhCpFitFile[2].data['I']
srhCorrV = srhCpFitFile[2].data['V']

srhDateTime = []
for srhT in srhTime[0]:
    srhHHMMSS = srhTime2HHMMSS(srhT)
    srhDateTime.append(datetime(int(dateList[0]),int(dateList[1]),int(dateList[2]),srhHHMMSS[0],srhHHMMSS[1],srhHHMMSS[2]))
    
fig = PL.figure(figsize=(12,10))
fig.suptitle(srhCpFitFile[0].header['DATE-OBS'])
pl = fig.subplots(nrows=2,ncols=1)

pl[0].plot(goesDatetime[smoothWindow:-smoothWindow],smooth(goesFile.variables[var_name[0]][:],smoothWindow)[smoothWindow:-smoothWindow],color='blue', label="{} {}".format(platform, var_name[0]))
pl[0].plot(goesDatetime[smoothWindow:-smoothWindow],smooth(goesFile.variables[var_name[1]][:],smoothWindow)[smoothWindow:-smoothWindow],color='red', label="{} {}".format(platform, var_name[1]))
pl[0].plot([srhDateTime[0],srhDateTime[-1]],[1e-7,1e-7],label = 'A class',color='gray')
pl[0].set_yscale("log")
pl[0].set_xlabel("Time [UT]")
pl[0].set_ylabel("X-Ray Flux [{}]".format(goesFile[var_name[0]].units))
pl[0].set_xlim(goesDatetime[0],goesDatetime[36000])
pl[0].legend(prop={'size': 12})
pl[0].grid()
    
tpl = pl[0].twinx()
tpl.set_ylabel('correlation');
tpl.plot(srhDateTime,srhCorrI.mean(axis=0), label='SRH')
tpl.legend(loc='lower right',prop={'size': 12})

pl[1].set_title('Stereo A')
pl[1].xaxis.set_major_formatter(PL.FuncFormatter(hm_format))
pl[1].xaxis.set_major_locator(PL.MultipleLocator(120))
pl[1].yaxis.set_major_formatter(PL.FuncFormatter(freq_format))
pl[1].set_xlim(0,600)
#pl[1].set_ylim(1,350)
pl[1].set_xlabel('Time UT')
pl[1].set_ylabel('Frequency MHz')
pl[1].imshow(stereoA['spectrum'].T,origin='lower',vmin=0, vmax=10., aspect=1/3)
