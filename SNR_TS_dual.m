Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\19-SNR\';
flienameRegx = 'TS_GS=*_count=*.bin';
filelist = dir([Directory,flienameRegx]);
[~,index] = sortrows({filelist.date}.');
filelist = filelist(index); 
clear index;
rate = 1e8/6.4;

noiFloorSampleRange = [112000, 116000];
avg = 100;
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
for j = 1:length(fileidx)
    psds = cell(1, avg);
    laserRefPs = zeros(1, avg)
    tempidx = fileidx{j}(1)-1;
    parfor i = fileidx{j}(1):fileidx{j}(avg)
        %gs1 = regexp(filelist(1).name, 'TS_GS=([0-9]+\.[0-9]*)+_count=[0-9]+\.csv', 'tokens')
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
        fclose(fileID);
    end
    avgpsd = psds{1};
    for i = 2:avg
        avgpsd = avgpsd + psds{i};
    end
    avgpsd = avgpsd/100;
    avgpsd = 10*log10(avgpsd*20);
    avgpsds{j} = avgpsd;
end

data = importdata([filelist(1).folder,'\',filelist(1).name]);
data = data.data;
VA = data(:,2);
[~, f] = periodogram(VA, rectwin(length(VA)), length(VA), 1/(data(2,1)-data(1,1)));

peakF = 122000;
span = 2000;
peakSearchRange = [peakF-span/2, peakF+span/2];
[~,PeakIndexs(1)] = min(abs(f-peakSearchRange(1)));
[~,PeakIndexs(2)] = min(abs(f-peakSearchRange(2)));

noiFloorSampleRange = [112000, 116000];
[~,Indexs(1)] = min(abs(f-noiFloorSampleRange(1)));
[~,Indexs(2)] = min(abs(f-noiFloorSampleRange(2)));

height = zeros(1, length(avgpsds));
noi = zeros(1, length(avgpsds));
for i = 1:length(avgpsds)
    p = avgpsds{i}(PeakIndexs(1):PeakIndexs(2));
    [~,loc] = findpeaks(p, 'MinPeakDistance',PeakIndexs(2)-PeakIndexs(1)-1);
    power = 10*log10(sum(10.^(p(loc-5:loc+5)/10)));
    height(i) = power;
    noi(i) = mean(avgpsds{i}(noiFloorSampleRange(1):noiFloorSampleRange(2)));
end
gapsizes = max(gapsizes)-gapsizes;