Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\49-SNR\';
experienceLOC = 26;

cPowerFlienameRegx = 'cPower_GS=*.bin';
cPowerFilelist = dir([Directory,cPowerFlienameRegx]);
[~,index] = sortrows({cPowerFilelist.date}.');
cPowerFilelist = cPowerFilelist(index); 
clear index;

powerFlienameRegx = 'power_GS=*.bin';
powerFilelist = dir([Directory,powerFlienameRegx]);
[~,index] = sortrows({powerFilelist.date}.');
powerFilelist = powerFilelist(index); 
clear index;

PSDoverP_FlienameRegx = 'PSDoverP_GS=*.bin';
PSDoverP_Filelist = dir([Directory,PSDoverP_FlienameRegx]);
[~,index] = sortrows({PSDoverP_Filelist.date}.');
PSDoverP_Filelist = PSDoverP_Filelist(index); 
clear index;

PSDoverPP_FlienameRegx = 'PSDoverPP_GS=*.bin';
PSDoverPP_Filelist = dir([Directory,PSDoverPP_FlienameRegx]);
[~,index] = sortrows({PSDoverPP_Filelist.date}.');
PSDoverPP_Filelist = PSDoverPP_Filelist(index); 
clear index;

PSDoverP_std_FlienameRegx = 'PSDoverP_std_GS=*.bin';
PSDoverP_std_Filelist = dir([Directory,PSDoverP_std_FlienameRegx]);
[~,index] = sortrows({PSDoverP_std_Filelist.date}.');
PSDoverP_std_Filelist = PSDoverP_std_Filelist(index); 
clear index;

PSDoverPP_std_FlienameRegx = 'PSDoverPP_std_GS=*.bin';
PSDoverPP_std_Filelist = dir([Directory,PSDoverPP_std_FlienameRegx]);
[~,index] = sortrows({PSDoverPP_std_Filelist.date}.');
PSDoverPP_std_Filelist = PSDoverPP_std_Filelist(index); 
clear index;

fileID = fopen([Directory, 'PSDfreq.bin']);
freq = fread(fileID, 'double');
fclose(fileID);

fileID = fopen([Directory, 'DN_new.bin']);
DN = fread(fileID, 'double');
fclose(fileID);

info = jsondecode(fileread([Directory,'info.json']));

peakF = 122037;
span = 500;
peakSearchRange = [peakF-span/2, peakF+span/2];
[~,PI1] = min(abs(freq-peakSearchRange(1)));
[~,PI2] = min(abs(freq-peakSearchRange(2)));

noiFloorSampleRange = [113000, 116000];
[~,NI1] = min(abs(freq-noiFloorSampleRange(1)));
[~,NI2] = min(abs(freq-noiFloorSampleRange(2)));

gapsizes = zeros(1, size(PSDoverP_Filelist,1));
steps = length(gapsizes);

for i =1:size(PSDoverP_Filelist, 1)
    temp = regexp(PSDoverP_Filelist(i).name, 'PSDoverP_GS=([0-9]+\.[0-9]*)\.bin', 'tokens');
    gapsizes(i) = str2num(temp{1}{1});
end

avg = 1000;
cPowers = zeros(1, steps);
for i =1:steps
    temp = regexp(cPowerFilelist(i).name, 'cPower_GS=([0-9]+\.[0-9]*)\.bin', 'tokens');
    fileID = fopen([cPowerFilelist(i).folder,'\',cPowerFilelist(i).name]);
    data = fread(fileID, 'double');
    cPowers(i) = mean(data);
    fclose(fileID);
end

powers = zeros(1, steps);
for i =1:steps
    temp = regexp(powerFilelist(i).name, 'power_GS=([0-9]+\.[0-9]*)\.bin', 'tokens');
    fileID = fopen([powerFilelist(i).folder,'\',powerFilelist(i).name]);
    data = fread(fileID, 'double');
    powers(i) = mean(data);
    fclose(fileID);
end


nois = zeros(1, length(gapsizes));
heights = zeros(1, length(gapsizes));
noi_stds = zeros(1, length(gapsizes));
height_stds = zeros(1, length(gapsizes));
SNRoverP = zeros(1, length(gapsizes));
for i = 1 : steps
    
    fileID = fopen([PSDoverP_Filelist(i).folder,'\',PSDoverP_Filelist(i).name]);
    data = fread(fileID, 'double');
    PSD = data *cPowers(i);
    fclose(fileID);
    PSD = PSD - DN;

    nois(i) = mean(PSD(NI1:NI2));
    
    fileID = fopen([PSDoverP_std_Filelist(i).folder,'\',PSDoverP_std_Filelist(i).name]);
    data = fread(fileID, 'double');
    std = data *cPowers(i);
    fclose(fileID);
    noi_stds(i) = sqrt(sum(std(NI1:NI2).^2)/(NI2-NI1+1)^2);
    
    fileID = fopen([PSDoverPP_Filelist(i).folder,'\',PSDoverPP_Filelist(i).name]);
    data = fread(fileID, 'double');
    PSD = data *cPowers(i)*cPowers(i);
    fclose(fileID);
    PSD = PSD - DN;
    PSD = PSD(PI1:PI2);
    [~,loc] = findpeaks(PSD, 'MinPeakDistance', PI2-PI1-1);
    heights(i) = sum(PSD(loc-5:loc+5));
    
    fileID = fopen([PSDoverPP_std_Filelist(i).folder,'\',PSDoverPP_std_Filelist(i).name]);
    data = fread(fileID, 'double');
    std = data *cPowers(i)*cPowers(i);
    std = std(PI1:PI2);
    fclose(fileID);
    height_stds(i) = sqrt(sum(std(loc-5:loc+5).^2));
end
SNR = heights./nois;
SNR_std = sqrt((1./nois).^2.*height_stds.^2 + (heights./nois.^2).^2.*noi_stds.^2);
gapsizes = max(gapsizes)-gapsizes;

% fit
Axfit = erfFit(info.MirrorA.positions.value, powers, info.MirrorA.angle.value, info.Mirror_tilt_angle.value, 6.3*1.68);
GSfit_A = 2*(Axfit.x0*1e3-info.MirrorA.positions.value)*sin(info.MirrorA.angle.value)*cos(info.Mirror_tilt_angle.value);
GSfit_A_overW = GSfit_A/Axfit.wx*1e-3;

%plot
[~,cut] = min(abs(ds-GSfit_A_overW(end)));
cut = cut+1;
noi_offset = -0.45;
height_offset = 2.1;

subplot(2,3,1);
height_logstds = 10*height_stds./(log(10)*heights);
errorbar(GSfit_A_overW, 10*log10(20*heights), height_logstds,'.');
hold on;
plot(ds(1:cut), 10*log10(20*peakHeightC(1:cut)));
title('Subplot 1: Peak area');
xlabel('Gapsize/beamWaist')
ylabel('Power(mdB)')

subplot(2,3,2);
noi_logstds = 10*noi_stds./(log(10)*nois);
errorbar(GSfit_A_overW, 10*log10(20*nois), noi_logstds, '.');
hold on;
plot(ds(1:cut), 10*log10(20*noiFloorC(1:cut)));
title('Subplot 2: Shot noise');
xlabel('Gapsize/beamWaist')
ylabel('Power(mdB)')

subplot(2,3,3);
errorbar(GSfit_A_overW, 10*log10(20*SNR), sqrt(noi_logstds.^2+height_logstds.^2),'.');
hold on;
plot(ds(1:cut), 10*log10(20*peakHeightC(1:cut)./noiFloorC(1:cut)));
title('Subplot 3: SNR');
xlabel('Gapsize/beamWaist')
ylabel('Power(mdB)')

subplot(2,3,4);
height_logstds = 10*height_stds./(log(10)*heights);
errorbar(GSfit_A_overW, 10*log10(20*heights)+height_offset, height_logstds,'.');
hold on;
plot(ds(1:cut), 10*log10(20*peakHeightC(1:cut)));
title(sprintf('Subplot 4: Peak area, offseted by %0.2f dB', height_offset));
xlabel('Gapsize/beamWaist')
ylabel('Power(mdB)')

subplot(2,3,5);
noi_logstds = 10*noi_stds./(log(10)*nois);
errorbar(GSfit_A_overW, 10*log10(20*nois)+noi_offset, noi_logstds, '.');
hold on;
plot(ds(1:cut), 10*log10(20*noiFloorC(1:cut)));
title(sprintf('Subplot 5: Shot Noise, offseted by %0.2f dB', noi_offset));
xlabel('Gapsize/beamWaist')
ylabel('Power(mdB)')

subplot(2,3,6);
errorbar(GSfit_A_overW, 10*log10(20*SNR)+height_offset-noi_offset, sqrt(noi_logstds.^2+height_logstds.^2),'.');
hold on;
plot(ds(1:cut), 10*log10(20*peakHeightC(1:cut)./noiFloorC(1:cut)));
title(sprintf('Subplot 6: SNR, offseted by %0.2f dB', height_offset-noi_offset));
xlabel('Gapsize/beamWaist')
ylabel('Power(mdB)')

set(gcf, 'Position',  [100, 100, 1350, 900])
