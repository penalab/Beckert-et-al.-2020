
itd = -200:20:200;
repet = 1;

bounds = nan(length(itd), 4);

for t = 1:length(itd)
    bounds(t, :) = define_peak_edges(itd(t), repet);
end

save model_bounds itd bounds

%%

function [bounds] = define_peak_edges(ITD, repet, dir_name)
if ~exist('repet', 'var')
    repet = 0;
end

r_txt = sprintf('%03d', repet);

% % Directory where SynchResults folder is found

if ~exist('dir_name', 'var')
    dir_name = 'f:\desktop\Code\Reproducibility\model_ICx\'; %% Can change this value %%
    %dir_name = 'K:\python\model_icx_codeonly';
end

eval(['cd ' dir_name '\SynchResults'])

% % Simulation parameters

%Stimulus ITD                           %% Can change this value %%
if ~exist('ITD', 'var')
    ITD = 0;
end

% % Load parameters

eval(['load ' dir_name '\SynchResults\data\SynchParms bestITD'])
% weights
Nn = length(bestITD);

W = zeros(Nn,Nn);
for n = 1:Nn
    W(n,:) = 12*exp(-.5*((bestITD - bestITD(n))/30).^2);
end

% % Load ICx spikes
eval(['load ' dir_name '\SynchResults\data\icx\SynchICx_ITD' num2str(ITD) '_' r_txt ' spikesICx Rate'])

% % plot mean rate

figure(1);clf;
plot(Rate,'o-')
ylabel('Firing rate (spikes/s)','fontsize',15)
xlabel('Best ITD (\mus)','fontsize',15)

%% Define main and side-peaks

check = 0;
while check == 0
    bounds = nan(1, 4);
    bounds(1) = input('main peak lower bound ');
    bounds(2) = input('main peak upper bound ');
    bounds(3) = input('side peak lower bound ');
    bounds(4) = input('side peak upper bound ');
    
    if bounds(2) - bounds(1) < bounds(4) - bounds(3)
        disp([num2str(abs(bounds(2) - bounds(1))), ' ' num2str(abs(bounds(4) - bounds(3)))])
        disp(num2str(bounds))
        disp('side peak is too large. please match bounds')
    else
        check = 1;
    end
    
end

end