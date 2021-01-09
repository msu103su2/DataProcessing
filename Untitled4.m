Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\20-SNR\';
flienameRegx = 'TS_GS=*_count=*.bin';
filelist = dir([Directory,flienameRegx]);
[~,index] = sortrows({filelist.date}.');
filelist = filelist(index); 
clear index;
rate = 1e8/6.4;

fileID = fopen([filelist(1).folder,'\',filelist(1).name]);
data = fread(fileID, [3,inf], 'double');
data = data';
VA = data(:,2);
[~, f] = periodogram(VA, rectwin(length(VA)), length(VA), rate);
fclose(fileID);

peakF = 122000;
span = 2000;
peakSearchRange = [peakF-span/2, peakF+span/2];
[~,PeakIndexs(1)] = min(abs(f-peakSearchRange(1)));
[~,PeakIndexs(2)] = min(abs(f-peakSearchRange(2)));
PI1 = PeakIndexs(1);
PI2 = PeakIndexs(2);

noiFloorSampleRange = [112000, 116000];
[~,Indexs(1)] = min(abs(f-noiFloorSampleRange(1)));
[~,Indexs(2)] = min(abs(f-noiFloorSampleRange(2)));
NI1 = Indexs(1);
NI2 = Indexs(2);

avg = 1000;
pbests = zeros(1, size(filelist,1));
gapsizes = [];
for i =1:size(filelist,1)
    temp = regexp(filelist(i).name, 'TS_GS=([0-9]+\.[0-9]*)_count=([0-9]+)\.bin', 'tokens');
    gapsize = str2num(temp{1}{1});
    I = find(gapsizes == gapsize);
    if isempty(I)
        gapsizes = [gapsizes, gapsize];
    end
end
fileidx = cell(1, length(gapsizes));
for i =1:size(filelist,1)
    temp = regexp(filelist(i).name, 'TS_GS=([0-9]+\.[0-9]*)_count=([0-9]+)\.bin', 'tokens');
    gapsize = str2num(temp{1}{1});
    I = find(gapsizes == gapsize);
    fileidx{I} = [fileidx{I}, i];
end

avgpsds = cell(1, length(gapsizes));
peaks = zeros(1, avg);
nois = zeros(1,avg);
for j = 1:length(fileidx)
    laserRefPs = zeros(1, avg);
    tempidx = fileidx{j}(1)-1;
    psds = cell(1, avg);
    parfor i = fileidx{j}(1):fileidx{j}(avg)
        temp = regexp(filelist(i).name, 'TS_GS=([0-9]+\.[0-9]*)_count=([0-9]+)\.bin', 'tokens');
        idx = str2num(temp{1}{2});
        gapsize = str2num(temp{1}{1});
        fileID = fopen([filelist(i).folder,'\',filelist(i).name]);
        data = fread(fileID, [3,inf], 'double');
        data = data';
        pbest = PDdualAC_balance(data);
        VA = data(:,1)*pbest;
        VB = data(:,2)*(2-pbest);
        laserRefPs(i) = mean(data(:,3));
        [newpsd, ~] = periodogram(VA-VB, rectwin(length(VA)), length(VA), rate);
        psds{i-tempidx} = newpsd;
        pbests(i) = pbest;
        
        p = newpsd(PI1:PI2);
        [~,loc] = findpeaks(p, 'MinPeakDistance',PI2-PI1-1);
        peaks(i) = 10*log10(sum(p(loc-5:loc+5))*20);
        nois(i) = 10*log10(mean(newpsd(NI1:NI2))*20);
        fclose(fileID);
    end
end


gapsizes = max(gapsizes)-gapsizes;