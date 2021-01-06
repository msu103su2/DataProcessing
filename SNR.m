Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\10-SNR\';
flienameRegx = 'j=*_SNR_D=*.csv';
filelist = dir([Directory,flienameRegx]);


noiFloorSampleRange = [112000, 116000];
peakF = 122000;
span = 2000;
peakSearchRange = [peakF-span/2, peakF+span/2];

data = importdata([filelist(1).folder,'\',filelist(1).name]);
data = data.data;
data(:,1) = data(:,1);
f = data(:,1);
p = data(:,2);
[~,Indexs(1)] = min(abs(f-noiFloorSampleRange(1)));
[~,Indexs(2)] = min(abs(f-noiFloorSampleRange(2)));
[~,PeakIndexs(1)] = min(abs(f-peakSearchRange(1)));
[~,PeakIndexs(2)] = min(abs(f-peakSearchRange(2)));
Indexs = Indexs + 1;
PeakIndexs = PeakIndexs +1;
noiOpts = detectImportOptions([filelist(1).folder,'\',filelist(1).name]);
noiOpts.DataLines = Indexs;
peakOpts = detectImportOptions([filelist(1).folder,'\',filelist(1).name]);
peakOpts.DataLines = PeakIndexs;

noiFloor = [];
opticalPower = [];
PeakHeight = [];
gapSize = [];

data = importdata([Directory,'DN.csv']);
data = data.data;
data(:,1) = data(:,1);
f = data(:,1);
p = data(:,2);
[~,a] = min(abs(f-noiFloorSampleRange(1)));
[~,b] = min(abs(f-noiFloorSampleRange(2)));
darkNoi = 10.^(p(a:b)/10)/20;

for i=1:size(filelist,1)
    if(strcmp(filelist(i).name,'BPDNR.csv'))
    else
        M = readmatrix([filelist(i).folder,'\',filelist(i).name],noiOpts);
        f = M(:,1);
        p = M(:,2);
        shotNoi = 10.^(p/10)/20 - darkNoi;
        
        %plot(f,p);
        %waitforbuttonpress;
        
        %{
        [pks,locs, w, ~] = findpeaks(p,'MinPeakProminence',5);
        if (size(pks,1) > 2) | (max(pks) > 10)
            continue
        end
        
        toExclude = abs(p-mean(p))>std(p);
        shotNoi(toExclude) = [];
        %}
        noiFloor = [noiFloor, mean(shotNoi)];
        %window = Indexs(2)-Indexs(1);
        %plot(f,p,'*',data(Indexs(1)-window:Indexs(2)+window,1),data(Indexs(1)-window:Indexs(2)+window,2));
        %waitforbuttonpress;
        M = readmatrix([filelist(i).folder,'\',filelist(i).name],peakOpts);
        f = M(:,1);
        p = M(:,2);

        [~,loc] = findpeaks(p,'MinPeakDistance',PeakIndexs(2)-PeakIndexs(1)-1);
        power = 10*log10(sum(10.^(p(loc-5:loc+5)/10)));
        PeakHeight = [PeakHeight, power];
        match = regexp(filelist(i).name,'SNR_D=([0-9]+.*[0-9]*).csv','tokens');
        gapSize = [gapSize, 13-str2double(match{1}{1})];
    end
end
noiFloor = 10*log10(noiFloor*20);
y = PeakHeight - noiFloor;

[gapSize,I] = sort(gapSize);
x = gapSize;
y = y(I);
PeakHeight = PeakHeight(I);
noiFloor = noiFloor(I);
%{
i = 1; j=2;
newX = [];
newPH = [];
PHstd = [];
newNF = [];
NFstd = [];
while i < size(x,2)
    while j <= size(x,2)
        if j == size(x,2)
            newX = [newX,x(i)];
            newPH = [newPH, mean(PeakHeight(i:j-1))];
            PHstd = [PHstd, std(PeakHeight(i:j-1))];
            newNF = [newNF, mean(noiFloor(i:j-1))];
            NFstd = [NFstd, std(noiFloor(i:j-1))];
            i = j;
            break
        else
            if x(j) == x(i)
                j = j+1;
            else
                newX = [newX,x(i)];
                newPH = [newPH, mean(PeakHeight(i:j-1))];
                PHstd = [PHstd, std(PeakHeight(i:j-1))];
                newNF = [newNF, mean(noiFloor(i:j-1))];
                NFstd = [NFstd, std(noiFloor(i:j-1))];
                i = j;
                j = j+1;
            end
        end
    end
end

x = newX;
PeakHeight = newPH;
noiFloor = newNF;
%}