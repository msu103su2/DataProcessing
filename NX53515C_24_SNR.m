Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\33-SNR\';
experienceLOC = 26;

PSDFlienameRegx = 'PSD_GS=*_count=*.bin';
PSDFilelist = dir([Directory,PSDFlienameRegx]);
[~,index] = sortrows({PSDFilelist.date}.');
PSDFilelist = PSDFilelist(index); 
clear index;

cPowerFlienameRegx = 'cPower_GS=*.bin';
cPowerFilelist = dir([Directory,cPowerFlienameRegx]);
[~,index] = sortrows({cPowerFilelist.date}.');
cPowerFilelist = cPowerFilelist(index); 
clear index;

fileID = fopen([Directory, 'PSDfreq.bin']);
freq = fread(fileID, 'double');
fclose(fileID);

fileID = fopen([Directory, 'DN_new.bin']);
DN = fread(fileID, 'double');
fclose(fileID);

noiFloorSampleRange = [112000, 116000];

peakF = 122037;
span = 500;
peakSearchRange = [peakF-span/2, peakF+span/2];
[~,PI1] = min(abs(freq-peakSearchRange(1)));
[~,PI2] = min(abs(freq-peakSearchRange(2)));

noiFloorSampleRange = [112000, 116000];
[~,NI1] = min(abs(freq-noiFloorSampleRange(1)));
[~,NI2] = min(abs(freq-noiFloorSampleRange(2)));

gapsizes = zeros(1, size(cPowerFilelist,1));

avg = 1000;
p = cell(1, length(gapsizes));
for i =1:size(cPowerFilelist,1)
    temp = regexp(cPowerFilelist(i).name, 'cPower_GS=([0-9]+\.[0-9]*)\.bin', 'tokens');
    gapsizes(i) = str2num(temp{1}{1});
    fileID = fopen([cPowerFilelist(i).folder,'\',cPowerFilelist(i).name]);
    data = fread(fileID, 'double');
    p{i} = data(1:min(avg, length(data)));
    fclose(fileID);
end

fileidx = cell(1, length(gapsizes));
for i =1:size(PSDFilelist, 1)
    temp = regexp(PSDFilelist(i).name, 'PSD_GS=([0-9]+\.[0-9]*)_count=([0-9]+)\.bin', 'tokens');
    gapsize = str2num(temp{1}{1});
    I = find(gapsizes == gapsize);
    fileidx{I} = [fileidx{I}, i];
end

avgpsds = cell(1, length(gapsizes));
nois = cell(1, length(gapsizes));
heights = cell(1, length(gapsizes));
SNRoverP = cell(1, length(gapsizes));
parfor j = 1 : length(fileidx)
    N_read = min(avg,length(fileidx{j}));
    SNRoverP{j} = zeros(1, N_read);
    nois{j} = zeros(1, N_read);
    heights{j} = zeros(1, N_read);
    for i = fileidx{j}(1):fileidx{j}(N_read)
        temp = regexp(PSDFilelist(i).name, 'PSD_GS=([0-9]+\.[0-9]*)_count=([0-9]+)\.bin', 'tokens');
        gapsize = str2num(temp{1}{1});
        fileID = fopen([PSDFilelist(i).folder,'\',PSDFilelist(i).name]);
        psd = fread(fileID, 'double');
        fclose(fileID);
        psd = psd - DN;
        
        noi = mean(psd(NI1:NI2));
        psd = psd(PI1:PI2);
        [~,loc] = findpeaks(psd, 'MinPeakDistance', PI2-PI1-1);
        if (abs(loc - experienceLOC) < 7)
            height = sum(psd(loc-5:loc+5));
            SNRoverP{j}(i - fileidx{j}(1) + 1) = height/(noi*p{j}(i - fileidx{j}(1) + 1));
            nois{j}(i - fileidx{j}(1) + 1) = noi;
            heights{j}(i - fileidx{j}(1) + 1) = height;
        else
            SNRoverP{j}(i - fileidx{j}(1) + 1) = nan;
            nois{j}(i - fileidx{j}(1) + 1) = nan;
            heights{j}(i - fileidx{j}(1) + 1) = nan;
        end    
    end
end

avgSNRoverP = zeros(1, length(gapsizes));
stdSNRoverP = zeros(1, length(gapsizes));
s = cell(1, length(gapsizes));
%{
heights_std_dB = zeros(1, length(gapsizes));
nois_std_dB = zeros(1, length(gapsizes));
heights_mean_dB = zeros(1, length(gapsizes));
nois_mean_dB = zeros(1, length(gapsizes));
%}
for j = 1 : length(p)
    p{j} = (p{j}(~isnan(nois{j})))';
    nois{j} = nois{j}(~isnan(nois{j}));
    heights{j} = heights{j}(~isnan(heights{j}));
    
    s{j} = (p{j}.*nois{j})./heights{j};
    avgSNRoverP(j) = mean(s{j});
    
    stdSNRoverP(j) = std(s{j})/sqrt(length(nois{j}));
    %{
    temp = mean(nois{j});
    nois_mean_dB(j) = 10*log10(20*temp);
    nois_std_dB(j) = 10*std(nois{j})/(temp*sqrt(length(nois{j}))*log(10));
    
    temp = mean(heights{j});
    heights_mean_dB(j) = 10*log10(20*temp);
    heights_std_dB(j) = 10*std(heights{j})/(temp*sqrt(length(heights{j}))*log(10));
    %}
end

gapsizes = max(gapsizes)-gapsizes;

for j = 1 : length(p)
    p{j} = p{j}(2:end);
    nois{j} = nois{j}(2:end);
    heights{j} = heights{j}(2:end);
end

avg_h = zeros(1,length(p));
std_h = zeros(1,length(p));
avg_n = zeros(1,length(p));
std_n = zeros(1,length(p));
for j = 1: length(p)
    avg_h(j) = mean(heights{j});
    std_h(j) = std(heights{j})/sqrt(length(p{j}));
    avg_n(j) = mean(nois{j});
    std_n(j) = std(nois{j})/sqrt(length(p{j}));
end