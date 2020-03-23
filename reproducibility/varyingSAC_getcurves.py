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
import pickle

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
    files = sorted(listdir(path))
    return TrialData, dur, name, files

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
    files = sorted(listdir(path))
    return TrialData, dur, name, files

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
    files = sorted(listdir(path))
    return TrialData, dur, name, files

# OT FF
def load_ot_ff():
    dur = 0.15
    name = 'ot'
    data = sio.loadmat('f:\desktop\WorkingFolder\distance\data\FreeField_OT_python.mat')
    data = data['data1']
    TrialData = []
    files = []
    for f, path in enumerate(data):
        tmp = data[f][0]
        for ix, tt in enumerate(tmp):
            tmp[ix] = [np.reshape(np.concatenate(t), [1, len(t)])  if np.size(t) > 0 else t for t in tt]
        TrialData.append(tmp)
        files.append(f)
    return TrialData, dur, name, files

# Get the integral of SAC peak out
def peak_integral(SAC, int_w):
    repro = np.sum(SAC[(len(SAC)/2 + 1) - int_w:(len(SAC)/2 + 1) + int_w])
    return repro    

#%% Load data 

# Here execute the desired data set
TrialData, dur, name, files = load_icls_itd()
TrialData, dur, name, files = load_icls_ild()
TrialData, dur, name, files = load_icls_ff()
TrialData, dur, name, files = load_ot_ff()

#%% Actual Work
##### make histograms for each data file

chdir('f:\desktop\WorkingFolder')

note = '_10ms'

# some variables
spikebin = 0.01 * br.second

hist = []
full_fr = []
x_spikes = []
y_spikes = []

for idx, t in enumerate(TrialData):
    hist.append([])
    for xx, da in enumerate(t):
        reps = len(da)
        y = [np.ones(np.size(s, 1)) * (tt + 1) for tt, s in enumerate(da) if s.size != 0]
        check_dim = map(lambda x: x == 0, map(np.size, da))
        check_dim = [c for c, chck in enumerate(check_dim) if chck]
        da = np.delete(da, check_dim)
        da = [s[0] for s in da]
        a = [np.shape(r) for r in da]
        full_fr.append(np.sum([b[0] for b in a]) / float(reps) / (dur - 0.05))
        try:
            x_spikes.append(np.concatenate(da))
        except:
            x_spikes.append(da)
        try:
            y_spikes.append(np.concatenate(y))  
        except:
            y_spikes.append(y)              
        try:
            hist[idx].append(make_PSTH(np.concatenate(da), dur * br.second, spikebin, reps, 0.1 * br.second))
        except:
            hist[idx].append(None)

cw = spikebin / br.second
del spikebin, idx, da, check_dim, c, chck, TrialData, r, a, b, y, tt, s

# % % run reproducibility

# variables
curve_time = 100

repro = DataFrame(columns = ['neuron', 'depvar', 'fr', 'curve'])

for neu, h in enumerate(hist):
    
    for dv, p in enumerate(h):
        try:
            t, SAC = autocor(p, reps, T=curve_time*br.ms, bin=cw * br.second)
            SAC = SAC * br.Hz ** 2
            sp = np.mean(p)
        except:
            strength_b = np.nan
            arr = np.arange(-(curve_time * br.ms), (curve_time * br.ms) , cw * br.second)
            arr = np.delete(arr, 0)
            arr[:] = np.nan
            SAC = arr
            sp = np.nan
            pass
        tmp = {'neuron': neu, 'depvar': dv, 'fr':sp, 'curve': SAC}
        repro = repro.append(tmp, ignore_index = True)
    
repro['fr'] = repro['fr'] / br.Hz
repro['x_spikes'] = x_spikes
repro['y_spikes'] = y_spikes
repro['full_fr'] = full_fr

#sio.savemat('repro_ff.mat', {'fr':fr, 'r_brette':r_brette, 'r_joris':r_joris})

curve = [c for c in repro['curve']]
sio.savemat('repro_' + name + '_curve' + note + '.mat', {'curve':curve, 'cw':cw, 'neuron':repro['neuron'].tolist(), 'depvar':repro['depvar'].tolist(), 'fr':repro['fr'].tolist(), 'full_fr':full_fr, 'reps':reps, 'dur':dur})

#df = repro.drop(columns = ['curve'])
#export_csv =df.to_csv (r'f:\desktop\WorkingFolder\repro_' + name + '_' + note + '.csv', index = None, header = True)

del cw, tmp, neu, dv, h, p, t, sp, SAC, dur, hist, reps, x_spikes, y_spikes, full_fr

files = pd.DataFrame({'neuron':range(len(files)), 'file':files})
repro = pd.merge(files, repro)
pickle.dump(repro, open('repro_' + name + note + '.pkl', 'wb'))