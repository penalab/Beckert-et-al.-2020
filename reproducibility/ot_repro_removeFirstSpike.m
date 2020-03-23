load('K:\WorkingFolder\distance\data\FreeField_OT_python.mat')

for i = 1:length(data1)
    for r = 1:size(data1{i}, 1)
        for c = 1:size(data1{i}, 2)
            if ~isempty(data1{i}{r,c})
                data1{i}{r,c}(1) = [];
            end
        end
    end
end

clear i c r

save('K:\WorkingFolder\distance\data\FreeField_OT_python_firstSpikeRemoved.mat')
