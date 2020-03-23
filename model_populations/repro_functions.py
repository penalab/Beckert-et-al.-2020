# -*- coding: utf-8 -*-
"""
Functions to perform reproducibility analysis on various data sets
"""

#%% importing and defining functions

#### import some stuff

import numpy as np
import brian2 as br
from pandas import DataFrame


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

def peak_integral(SAC, int_w):
    repro = np.sum(SAC[(len(SAC)/2 + 1) - int_w:(len(SAC)/2 + 1) + int_w])
    return repro    

def step1_makehist_gettimes(TrialData, cw, dur):
    cw = cw * br.second
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
                hist[idx].append(make_PSTH(np.concatenate(da), dur * br.second, cw, reps, 0.1 * br.second))
            except:
                hist[idx].append(None)
    return hist, full_fr, x_spikes, y_spikes, reps

def step2_calcSAC(hist, cw, reps, curve_time = 100):
    repro = DataFrame(columns = ['neuron', 'depvar', 'fr', 'curve'])
    for neu, h in enumerate(hist):
        for dv, p in enumerate(h):
            try:
                t, SAC = autocor(p, reps, T=curve_time*br.ms, bin=cw * br.second)
                SAC = SAC * br.Hz ** 2
                sp = np.mean(p)
            except:
                arr = np.arange(-(curve_time * br.ms), (curve_time * br.ms) , cw * br.second)
                arr = np.delete(arr, 0)
                arr[:] = np.nan
                SAC = arr
                sp = np.nan
                pass
            tmp = {'neuron': neu, 'depvar': dv, 'fr':sp, 'curve': SAC}
            repro = repro.append(tmp, ignore_index = True)
    return repro

def step2_calcSAC_stripped(hist, cw, reps, curve_time = 100):
    repro = DataFrame(columns = ['neuron', 'depvar', 'fr', 'curve'])
    for neu, p in enumerate(hist):
        try:
            t, SAC = autocor(p, reps, T=curve_time*br.ms, bin=cw * br.second)
            SAC = SAC * br.Hz ** 2
            sp = np.mean(p)
        except:
            arr = np.arange(-(curve_time * br.ms), (curve_time * br.ms) , cw * br.second)
            arr = np.delete(arr, 0)
            arr[:] = np.nan
            SAC = arr
            sp = np.nan
            pass
        tmp = {'neuron': neu, 'depvar': 0, 'fr':sp, 'curve': SAC}
        repro = repro.append(tmp, ignore_index = True)
    return repro

def best_fr(hist):
    hist_strip = []
    for neu, h in enumerate(hist):
        fr = [np.mean(t) if t is not None else np.nan for t in h]
        try:
            ix = np.nanargmax(fr)
            hist_strip.append(hist[neu][ix])
        except:
            print('error with neuron %s' %neu)    
    return hist_strip
    
    
    