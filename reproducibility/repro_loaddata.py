# -*- coding: utf-8 -*-
"""
Some automated functions to extract the raw data to be analyzed
"""

from os import chdir, listdir
import scipy.io as sio
import numpy as np

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
    chdir('K:\Data\ICls\Data\RawData\Dichotic\ITD')
    path= 'K:\Data\ICls\Data\RawData\Dichotic\ITD'
    
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
    chdir('K:\Data\ICls\Data\RawData\Dichotic\ILD')
    path= 'K:\Data\ICls\Data\RawData\Dichotic\ILD'
    
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
    chdir('K:\Data\ICls\Data\RawData\FF')
    path = 'K:\Data\ICls\Data\RawData\FF'
    
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
    data = sio.loadmat('K:\WorkingFolder\distance\data\FreeField_OT_python.mat')
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