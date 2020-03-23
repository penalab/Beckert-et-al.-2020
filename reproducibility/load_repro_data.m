FR = cell(6, 1);
REP = cell(6, 1);

for i = 1:4
    switch i
        case 1
            load('K:\WorkingFolder\working_reprod\frompython\repro_itd_curve_100micro.mat')
            mw_time = 0.0005;
        case 2
            load('K:\WorkingFolder\working_reprod\frompython\repro_ild_curve_100micro.mat')
            mw_time = 0.0005;
        case 3
            load('K:\WorkingFolder\working_reprod\frompython\repro_az_curve_100micro.mat')
            mw_time = 0.0005;
        case 4
            load('K:\WorkingFolder\working_reprod\frompython\repro_el_curve_100micro.mat')
            mw_time = 0.0005;
    end
    
    center = round(size(curve, 2) / 2);
    win = round(mw_time / cw);
    winsac = curve(:, center - win + 1:center + win - 1);
    f = repmat(fr' .^2, 1, size(winsac, 2));
    repro = sum(winsac - f, 2) * cw ./ fr';
    
    neu = unique(neuron);
    FR{i} = cell(length(neu), 1);
    REP{i} = cell(length(neu), 1);
    
for n = 1:length(neu)
    idx = neuron == neu(n);
    REP{i}{n} = repro(idx)';
    FR{i}{n} = fr(idx);
end  
    
end

clear i center win winsac f repro neu n idx mw_time curve cw neuron fr depvar ex_ix e full_fr

d = load('K:\WorkingFolder\working_reprod\frompython\repro_otaz_curve_1ms.mat');
FR{5} = d.FR;
REP{5} = d.REP;
d = load('K:\WorkingFolder\working_reprod\frompython\repro_otel_curve_1ms.mat');
FR{6} = d.FR;
REP{6} = d.REP;

clear d