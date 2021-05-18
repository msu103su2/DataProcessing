Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\72-CM_test\';

FlienameRegx = 'f=*.bin';
Filelist = dir([Directory,FlienameRegx]);
[~,index] = sortrows({Filelist.date}.');
Filelist = Filelist(index); 
clear index;

fileID = fopen([Directory, 'PSDfreq.bin']);
freq = fread(fileID, 'double');
fclose(fileID);

N = length(Filelist);
sumSpan = 50;

fs = zeros(1, N);
Is = zeros(1, N);
Powers = zeros(1, N);
for i = 1 : N
    temp = regexp(Filelist(i).name, 'f=([0-9]+\.[0-9]*)\.bin', 'tokens');
    fs(i) = str2num(temp{1}{1});
    [~,Is(i)] = min(abs(freq-fs(i)));
    
    fileID = fopen([Filelist(i).folder,'\',Filelist(i).name]);
    data = fread(fileID, 'double');
    fclose(fileID);
    %plot(freq(Is(i)-1000:Is(i)+1000),10*log10(20*data(Is(i)-1000:Is(i)+1000)));
    %waitforbuttonpress;
    
    Powers(i) = sum(data(Is(i)-sumSpan:Is(i)+sumSpan));
end

plot(fs, Powers);

