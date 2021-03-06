Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\19-SNR\';
flienameRegx = 'TS_GS=*_count=*.csv';
filelist = dir([Directory,flienameRegx]);
[~,index] = sortrows({filelist.date}.');
filelist = filelist(index); 
clear index;
avg = 1000;
gapsizes = [];
for i =1:size(filelist,1)
    temp = regexp(filelist(i).name, 'TS_GS=([0-9]+\.[0-9]*)_count=([0-9]+)\.csv', 'tokens');
    gapsize = str2num(temp{1}{1});
    I = find(gapsizes == gapsize);
    if isempty(I)
        gapsizes = [gapsizes, gapsize];
    end
end
fileidx = cell(1, length(gapsizes));
for i =1:size(filelist,1)
    temp = regexp(filelist(i).name, 'TS_GS=([0-9]+\.[0-9]*)_count=([0-9]+)\.csv', 'tokens');
    gapsize = str2num(temp{1}{1});
    I = find(gapsizes == gapsize);
    fileidx{I} = [fileidx{I}, i];
end

powers = zeros(1, length(gapsizes));
for j = 1:length(fileidx)
    power = zeros(1, avg);
    for i = fileidx{j}(1):fileidx{j}(avg)
        temp = regexp(filelist(i).name, 'TS_GS=([0-9]+\.[0-9]*)_count=([0-9]+)\.csv', 'tokens');
        data = importdata([filelist(i).folder,'\',filelist(i).name]);
        data = data.data;
        VC = data(:,3)
        power(i) = mean(VC);
    end
    powers(j) = mean(power);
end