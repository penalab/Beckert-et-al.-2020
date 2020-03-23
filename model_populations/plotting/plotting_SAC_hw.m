% c = [];

for i = 1:5
    switch i
        case 1
            load('K:\WorkingFolder\working_reprod\frompython\repro_itd_curve_100micro.mat', 'curve')
        case 2
            load('K:\WorkingFolder\working_reprod\frompython\repro_ild_curve_100micro.mat', 'curve')
        case 3
            load('K:\WorkingFolder\working_reprod\frompython\repro_az_curve_100micro.mat', 'curve')
        case 4
            load('K:\WorkingFolder\working_reprod\frompython\repro_el_curve_100micro.mat', 'curve')
        case 5
            load('K:\WorkingFolder\working_reprod\frompython\repro_ot_curve_100micro.mat', 'curve')
    end
    

   c = [c; curve];
    
end

c = norm01(c);

clear curve i

xx = -100:0.1:100;
xx = xx(2:end-1);

plot(xx, nanmean(c)/max(nanmean(c)), 'r');

% load('K:\WorkingFolder\working_reprod\frompython\repro_ot_curve_100micro.mat', 'curve');
% c_ot = curve;
% hold on 
% plot(xx, nanmean(c_ot)/max(nanmean(c_ot)), 'b');