function power_position_check()
    %the function is used to check wether the measured average laser power
    %from one of the split mirror is in good agreement with its measured position
    MaxV_total = 17.95; %max volt measured from both PDs
    
    Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\53-SNR\';
    info = jsondecode(fileread([Directory,'info.json']));
    
    powerFlienameRegx = 'power_GS=*.bin';
    powerFilelist = dir([Directory,powerFlienameRegx]);
    [~,index] = sortrows({powerFilelist.date}.');
    powerFilelist = powerFilelist(index); 
    clear index;
    
    cPowerFlienameRegx = 'cPower_GS=*.bin';
    cPowerFilelist = dir([Directory,cPowerFlienameRegx]);
    [~,index] = sortrows({cPowerFilelist.date}.');
    cPowerFilelist = cPowerFilelist(index); 
    clear index;
    
    steps = length(powerFilelist);
    n_powers = zeros(1, steps);
    for i =1:steps
        fileID = fopen([powerFilelist(i).folder,'\',powerFilelist(i).name]);
        data = fread(fileID, 'double');
        fclose(fileID);

        fileID = fopen([cPowerFilelist(i).folder,'\',cPowerFilelist(i).name]);
        data2 = fread(fileID, 'double');
        fclose(fileID);

        n_powers(i) = mean(data./data2)*mean(data2); %as laser fluctuate in time, average by normalizing to it
    end
    
    [param, ~] = fitOptions('B');
    x = (param(3)*1e3-info.MirrorB.positions.value)*sin(info.MirrorB.angle.value)*cos(info.Mirror_tilt_angle.value);
    x_overW = x/param(1)*1e-3;
    relativePower = (n_powers)/(MaxV_total);
    
    relativePower_T = arrayfun(@(x) func(x), -x_overW);
    
    plot(x_overW, relativePower_T, x_overW, relativePower, '*');
end

function re = func(x)
    re = 1/2*(1+erf(x*sqrt(2)));
end

function [param, param_std] = fitOptions(Mirror)
    Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\55-Cali\';
    CaliFlienameRegx = 'Cali*.json';
    CaliFilelist = dir([Directory,CaliFlienameRegx]);
    
    x0 = [];
    w0 = [];
    for i = 1 : length(CaliFilelist)
        fitInfo = jsondecode(fileread([Directory,CaliFilelist(i).name]));
        if fitInfo.CaliMirror == Mirror
            x0 = [x0, fitInfo.fitResult(1)];
            w0 = [w0, fitInfo.fitResult(end)];
        end
    end
    param(3) = mean(x0);
    param_std(3) = std(x0);
    param(1) = mean(w0);
    param_std(1) = std(w0);
end