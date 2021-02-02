Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\20-SNR\';
flienameRegx = 'TS_GS=*_count=*.bin';
filelist = dir([Directory,flienameRegx]);
[~,index] = sortrows({filelist.date}.');
filelist = filelist(index); 
clear index;
rate = 1e8/6.4;

checkrange = [110e3, 130e3];
psds = cell(1, 100);
parfor i = 1:100
    fileID = fopen([filelist(i).folder,'\',filelist(i).name]);
    data = fread(fileID, [3,inf], 'double');
    data = data';
    VA = data(:,1)*pypbests(i);
    VB = data(:,2)*(2-pypbests(i));
    [newpsd, f] = periodogram(VA-VB, rectwin(length(VA)), length(VA), rate);
    psds{i} = newpsd;
    fclose(fileID);
end

pypsd = zeros(size(psds{1}));
for i = 1:100
    pypsd = pypsd+psds{i};
end
pypsd = 10*log10(20*pypsd/100);