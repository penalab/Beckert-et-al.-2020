function bsRaster(spikes,bins,stim)
% bsRaster                  Plot a raster.
%
% BSRASTER(SPIKES,BINS,STIM) will plot a raster of the
% spikes in SPIKES between the bounds specified by BINS (of
% the form [lower upper]) relative to the time STIM (all
% times in ms).  STIM can be either a scalar (which causes
% all trials to be aligned to the same time) or a vector of
% the same length as the number of trials.  SPIKES is a cell
% matrix, rows being different cells, columns being
% different trials, and each cell being an array of spike
% times (consistent with the output of BEDSgetSpikes).
%
% By default, BINS is [0 500] and STIM is 0.
%
% For unknown reasons, attempting to plot the rasters of
% more than 50 trials at the same time will cause Matlab
% 6.0 to hang.  Matlab 5.3 does not have this problem.
%
% See also:  BEDS, bsPsth.

% RCS:$Id: bsRaster.m,v 1.4 2001/10/13 17:37:57 bjorn Exp $
% Written by GBC, 21 May 2001.
% Copyright 1997, California Institute of Technology.
%	    written  by Maneesh Sahani (maneesh@caltech.edu)
%	    modified by John Pezaris   (pz@caltech.edu)

% default the arguments
color = num2cell(get(gca, 'colororder'), 2);

if (nargin < 2) bins = [0 500]; end
if (nargin < 3) stim = 0;       end

% Make sure stim is the right size.
if length(stim) == 1
  stim = stim * ones(size(spikes(:,1)));
end

if length(stim) > 0
  numtrials = length(spikes(:,1));
  numcells = length(spikes(1,:));

  for c = 1:numcells
    cell_color = color{2+mod(c-1, length(color))};
    
    for t = 1:numtrials
      evt = spikes{t,c} - stim(t);
      if length(evt)
	valevt = find((evt > bins(1)) & (evt < bins(2)));
      end
      if length(evt)
	set(gca,'ColorOrder',cell_color,'NextPlot','add'); 
	h = plot(evt(valevt), (t+(c-1)*numtrials)*ones(size(evt(valevt))), 'k.');
	set(h, 'MarkerSize', 5.5);   % reduce dot size
      end
    end
  end    
  set(gca, 'ylim', [0 numcells * numtrials+1]);
end


  




