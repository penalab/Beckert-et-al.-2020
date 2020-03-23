import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import brian2 as br
from os import chdir, listdir

chdir('Z:\Michael Beckert\data\ICls\Data\RawData\Dichotic\ILD')
path = 'Z:\Michael Beckert\data\ICls\Data\RawData\Dichotic\ILD'
files = sorted(listdir(path))

chdir('f:\desktop\WorkingFolder')

repro = pd.read_pickle('./python_repro/repro_ild_base_50micro.pkl')

#files = pd.DataFrame({'neuron':range(len(files)), 'file':files})
#repro = pd.merge(files, repro)

del path, files
#%%

cw = 0.00005 * br.second
mw = 0.05 * br.second
times = np.arange(-mw + cw, mw - 0.00000001 * br.second, cw)

curves = [(c - (f * br.Hz) ** 2) * cw / (f * br.Hz) for c, f in zip(repro['curve'], repro['fr'])]
reli = pd.Series([c[len(c) / 2 + 1] for c in curves])

for n in np.unique(repro['neuron']):
#n = 63

    ix = np.where(repro['neuron'] == n)
    axes = plt.subplots(squeeze = False)
    
    for count, i in enumerate(ix[0]):
        
        c = (repro['curve'][i] - (repro['fr'][i] * br.Hz) ** 2) * cw / (repro['fr'][i] * br.Hz)
        
        axes = plt.subplot2grid((np.size(ix, 1) + 1, 2), (count, 0))
        axes.plot(times, c)
        axes.set(ylabel = 'normalized coincidences')
        if count == len(ix[0]) - 1:
            axes.set(xlabel = 'lag (s)')
        
        axes = plt.subplot2grid((np.size(ix, 1) + 1, 2), (count, 1))
            
        axes.scatter(repro['x_spikes'][i], repro['y_spikes'][i], 0.5)
        axes.set(ylabel = 'trials')
        if count == len(ix[0]) - 1:
            axes.set(xlabel = 'time (s)')
             
        plt.show()
    
    r = reli[np.concatenate(ix)]
    f = repro['fr'][np.concatenate(ix)]
    
    axes = plt.subplot2grid((np.size(ix, 1) + 1, 2), (np.size(ix, 1), 0))
        
    axes.scatter(f, r)
    axes.set_xlim([np.min(f), np.max(f)])
    axes.set_ylim([np.min(r), np.max(r)])
    
    plt.suptitle(str(n) + '\n\n' + repro['file'][ix[0][0]])

del n, ix, axes, count, i, cw, mw, times, f, r, curves, reli, c