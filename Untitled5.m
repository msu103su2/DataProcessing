Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\53-SNR\';

PSDFlienameRegx = 'PSD_GS=*_count=*.bin';
PSDFilelist = dir([Directory,PSDFlienameRegx]);
[~,index] = sortrows({PSDFilelist.date}.');
PSDFilelist = PSDFilelist(index); 
%PSDFilelist = PSDFilelist(1:100); 
clear index;

fileID = fopen([Directory, 'PSDfreq.bin']);
freq = fread(fileID, 'double');
fclose(fileID);


peakF = 122037;
span = 500;
peakSearchRange = [peakF-span/2, peakF+span/2];
[~,PI1] = min(abs(freq-peakSearchRange(1)));
[~,PI2] = min(abs(freq-peakSearchRange(2)));

experienceLOC = 26;

areas = zeros(1, length(PSDFilelist));
for i = 1: length(PSDFilelist)
    fileID = fopen([PSDFilelist(i).folder,'\',PSDFilelist(i).name]);
    psd = fread(fileID, 'double');
    fclose(fileID);
    
    psd = psd(PI1:PI2);
    psd = psd(experienceLOC - 5:experienceLOC + 5);
    areas(i) = sum(psd);    
end
std1 = std(areas)/sqrt(length(PSDFilelist));

bins = zeros(length(PSDFilelist),11);
for i = 1: length(PSDFilelist)
    fileID = fopen([PSDFilelist(i).folder,'\',PSDFilelist(i).name]);
    psd = fread(fileID, 'double');
    fclose(fileID);
    
    psd = psd(PI1:PI2);
    bins(i,:) = psd(experienceLOC - 5:experienceLOC + 5);   
end
binstds = zeros(1,11);
for i = 1:11
    binstds(i) = std(bins(:,i));
end
std2 = sqrt(sum(binstds.^2))/sqrt(length(PSDFilelist));

PSD = zeros(size(freq));
PSD2 = zeros(size(freq));
for i = 1: length(PSDFilelist)
    fileID = fopen([PSDFilelist(i).folder,'\',PSDFilelist(i).name]);
    new = fread(fileID, 'double');
    PSD = PSD + new;
    PSD2 = PSD2 + new.^2;
    fclose(fileID);
end
N = length(PSDFilelist);
PSD = PSD/N;
std = sqrt((PSD2/N-PSD.^2))/sqrt(N);