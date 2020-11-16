f1 = figure;
folder = 'Z:\data\optical lever project\doubleRec\100_test\';
flienameRegx = '*.csv';
filelist = dir([folder,flienameRegx]);
locs = zeros(size(filelist,1),1);
time = repmat(datetime('18-Oct-2020 20:25:43'), size(filelist,1),1);
for i = 1: size(filelist,1)    
    table = importfile_Tshift([folder,filelist(i).name]);
    [~,locs(i)] = findpeaks(table.PSDdB2Hz,table.FrequenciesHz,'MinPeakDistance', table.FrequenciesHz(end)-table.FrequenciesHz(1)-1);
    time(i) = datetime(filelist(i).date);
end