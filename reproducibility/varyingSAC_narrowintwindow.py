# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 10:48:39 2018

@author: penalab
"""

#%% importing and defining functions

#### import some stuff

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from os import chdir, listdir
import scipy.io as sio
import numpy as np
import brian2 as br
from pandas import Series, DataFrame
import pandas as pd

# % % define some useful functions

def flatten_array(data):
    d = data[0][0].flatten()
    for i in np.arange(0, len(d)):
        if d[i] is not np.empty:
            d[i] = d[i][0]
        else:
            d[i] = np.empty(1)   
    out = np.concatenate(([i for i in d]))
    return out

def autocor(PSTH,N=None,T=20*br.ms,bin=None):
    if bin is None:
        bin = br.defaultclock.dt
    p = int(T/ bin)
    SAC = br.zeros(p)    
    if N is None:
        SAC[0] = br.mean(PSTH * PSTH)
    else: # correction to exclude self-coincidences
        PSTHnoself = br.clip(PSTH - 1. / (bin * N), 0*br.kHz, br.Inf*br.kHz)
        SAC[0] = br.mean(PSTH * PSTHnoself) * N / (N - 1.)
    SAC[1:] = [br.mean(PSTH[:-i] * PSTH[i:]) for i in range(1,p)]
    SAC = br.hstack((SAC[::-1], SAC[1:]))
    out = (br.arange(len(SAC)) - len(SAC) / 2) * bin    
    return out, SAC

def make_PSTH(spikes_flat, dur, cw, reps, start):
    PSTH, bin_edges = np.histogram(spikes_flat, bins = np.arange(start, dur + 0.05 * br.second + cw, cw))
    PSTH = PSTH / (cw * reps)
    return PSTH

def jorisnorm(sac, N, m, dur, cw):
    nf = N * (N - 1) * m**2 * dur * cw * br.hertz ** 2
    sacnorm = sac / nf
    strength = np.mean(sacnorm)
    return sacnorm, strength

def brettenorm(sac, m, cw):
    sacnorm = (sac - m ** 2) * cw / m
    strength = sum(sac - m ** 2) * cw / m
    return sacnorm, strength

def flatten(l):
  out = []
  for item in l:
    if isinstance(item, (list, tuple)):
      out.extend(flatten(item))
    else:
      out.append(item)
  return out

def convert_ms2sec( mat ):
    for r, c in np.ndindex(mat.shape): 
        mat[r][c] = mat[r][c] / 1000
    return mat

# % % Load in data from MatLab files
  # using functions to ease the amount of necessary commenting
    
# ICls ITD
def load_icls_itd():
    dur = 0.25
    name = 'itd'
    chdir('Z:\Michael Beckert\data\ICls\Data\RawData\Dichotic\ITD')
    path= 'Z:\Michael Beckert\data\ICls\Data\RawData\Dichotic\ITD'
    
    TrialData = []
    for f in sorted(listdir(path)):
        data = sio.loadmat(f)
        curvedata = data['curvedata']
        TrialData.append(curvedata['spike_times'][0][0])
        del data, curvedata
    TrialData = map(convert_ms2sec, TrialData)
    return TrialData, dur, name

# ICls ILD
def load_icls_ild():    
    dur = 0.25
    name = 'ild'
    chdir('Z:\Michael Beckert\data\ICls\Data\RawData\Dichotic\ILD')
    path= 'Z:\Michael Beckert\data\ICls\Data\RawData\Dichotic\ILD'
    
    TrialData = []
    for f in sorted(listdir(path)):
        data = sio.loadmat(f)
        curvedata = data['curvedata']
        TrialData.append(curvedata['spike_times'][0][0])
        del data, curvedata
    TrialData = map(convert_ms2sec, TrialData)
    return TrialData, dur, name

# ICls FF
def load_icls_ff():
    dur = 0.45
    name = 'ff'
    chdir('Z:\Michael Beckert\data\ICls\Data\RawData\FF')
    path = 'Z:\Michael Beckert\data\ICls\Data\RawData\FF'
    
    TrialData = []
    for f in sorted(listdir(path)):
        data = sio.loadmat(f)
        TrialData.append(data['TrialData'])
        del data
    TrialData = map(convert_ms2sec, TrialData)
    return TrialData, dur, name

# OT FF
def load_ot_ff():
    dur = 0.15
    name = 'ot'
    data = sio.loadmat('f:\desktop\WorkingFolder\distance\data\FreeField_OT_python.mat')
    data = data['data1']
    TrialData = []
    for f, path in enumerate(data):
        TrialData.append(data[f][0])
    return TrialData, dur, name

# Get the integral of SAC peak out
def peak_integral(SAC, int_w):
    repro = np.sum(SAC[(len(SAC)/2 + 1) - int_w:(len(SAC)/2 + 1) + int_w])
    return repro    

#%% Load data 

# Here execute the desired data set
TrialData, dur, name = load_icls_itd()
TrialData, dur, name = load_icls_ild()
TrialData, dur, name = load_icls_ff()
TrialData, dur, name = load_ot_ff()

#%% Actual Work
##### make histograms for each data file

chdir('f:\desktop\WorkingFolder')

note = 'base'

# some variables
spikebin = 0.0001 * br.second

hist = []

for idx, t in enumerate(TrialData):
    hist.append([])
    for da in t:
        reps = len(da)
        check_dim = map(lambda x: x != 1, map(len, da))
        check_dim = [c for c, chck in enumerate(check_dim) if chck]
        da = np.delete(da, check_dim)
        try:
            hist[idx].append(make_PSTH(np.concatenate(da, 1), dur * br.second, spikebin, reps, 0.1 * br.second))
        except:
            hist[idx].append(None)

del spikebin, idx, t, da, check_dim, c, chck, TrialData

# % % run reproducibility

# variables
cw = 0.0001
mw_time = 10
curve_time = 100
maxWidth = int(round((mw_time * br.ms) / (cw * br.second)))
bins = np.arange(0, mw_time, cw)

r_brette = DataFrame(columns = ['neuron', 'depvar', 'fr', 'curve', 'strength'])
r_joris = DataFrame(columns = ['neuron', 'depvar', 'curve', 'strength'])
fr = DataFrame(columns = ['neuron', 'depvar', 'fr'])

for neu, h in enumerate(hist):
    
    for dv, p in enumerate(h):
        try:
            t, SAC = autocor(p, reps, T=curve_time*br.ms, bin=cw * br.second)
            # t, SAC = autocor(p, reps, T=30*br.ms)
            SAC = SAC * br.Hz ** 2
            sacnorm_j, tmp = jorisnorm(SAC, reps, np.mean(p), dur, cw)
            sacnorm_b, tmp = brettenorm(SAC, np.mean(p), cw * br.second)
            SAC = SAC[(len(SAC)/2 + 1) - maxWidth:(len(SAC)/2 + 1) + maxWidth - 1]
            tmp, strength_j = jorisnorm(SAC, reps, np.mean(p), dur, cw)
            tmp, strength_b = brettenorm(SAC, np.mean(p), cw * br.second)
            sp = np.mean(p)
            tmp = {'neuron': neu, 'depvar': dv, 'fr': sp}
            fr = fr.append(tmp, ignore_index = True)
        except:
            strength_j = np.nan
            strength_b = np.nan
            arr = np.arange(-(curve_time * br.ms), (curve_time * br.ms) , cw * br.second)
            arr = np.delete(arr, 0)
            arr[:] = np.nan
            sacnorm_j = arr
            sacnorm_b = arr
            sp = np.nan
            tmp = {'neuron': neu, 'depvar': dv, 'fr': np.nan}
            fr.append(tmp, ignore_index = True)
            pass
        tmp = {'neuron': neu, 'depvar': dv, 'fr':sp, 'curve': sacnorm_b, 'strength': strength_b}
        r_brette = r_brette.append(tmp, ignore_index = True)
        tmp = {'neuron': neu, 'depvar': dv, 'curve': sacnorm_j, 'strength': strength_j}
        r_joris = r_joris.append(tmp, ignore_index = True)
    
del cw, tmp, maxWidth, bins, neu, dv, h, p, t, sp, SAC, strength_j, strength_b, dur, hist, reps, sacnorm_b, sacnorm_j

r_brette['strength'] = r_brette['strength'] / br.Hz
r_brette['fr'] = r_brette['fr'] / br.Hz

#sio.savemat('repro_ff.mat', {'fr':fr, 'r_brette':r_brette, 'r_joris':r_joris})

#curve = [c for c in r_brette['curve']]
#sio.savemat('repro_' + name + '_curve' + note + '.mat', {'curve':curve})

df = r_brette.drop(columns = ['curve'])
export_csv =df.to_csv (r'f:\desktop\WorkingFolder\repro_' + name + '_' + str(mw_time) + note + '.csv', index = None, header = True)



#fr = flatten(fr)
#r_brette = flatten(r_brette)
#r_joris = flatten(r_joris)
#
## % % Plot summary data
#
#red_patch = mpatches.Patch(color = 'red', label = 'Brette')
#blue_patch = mpatches.Patch(color = 'blue', label = 'Joris')
#
#fig, ax1 = plt.subplots()
#ax2 = ax1.twinx()
#ax1.scatter(fr, r_brette, c = 'r')
#ax2.scatter(fr, r_joris, c = 'b')
#ax1.set_xlabel('firing rate')
#ax1.set_ylabel('brette', color='r')
#ax2.set_ylabel('joris', color='b')
#plt.legend(handles = [red_patch, blue_patch])
#plt.title('Reproducbility vs Firing Rate for both types of normalization')
#
#plt.show