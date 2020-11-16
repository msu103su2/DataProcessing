folder = 'Z:/data/optical lever project/NORCADA_NX53515C/02-balancing/';
flienameRegx = 'TBPump20*.csv'
global samplesPerFile;
samplesPerFile = 801;
filelist = dir([folder,flienameRegx]);
dataTable = table('Size',[samplesPerFile*size(filelist(1).name,1) 2],'VariableNames', {'FrequenciesHz', 'PSDdB2Hz'},'VariableTypes', {'double', 'double'});
allModeTable = table('Size',[0 3],'VariableNames', {'relativeHeights', 'Heights', 'pksLocation'},'VariableTypes', {'double', 'double', 'double'});
for i = 1 : size(filelist,1)
    newtable = importfile1([folder,filelist(i).name]);
    dataTable((i-1)*samplesPerFile + 1 : i * samplesPerFile ,{'FrequenciesHz', 'PSDdB2Hz'}) = newtable(:,{'FrequenciesHz', 'PSDdB2Hz'});
    modeTable = FindPeaks(newtable.FrequenciesHz, newtable.PSDdB2Hz);
    allModeTable = [allModeTable; modeTable];
end
figure('Name',flienameRegx);
plot(dataTable.FrequenciesHz, dataTable.PSDdB2Hz, allModeTable.pksLocation, allModeTable.Heights, '*');

function modeTable = FindPeaks(X,Y)
    Mean = mean(Y);
    Std = std(Y);
    [Heights,pksLocation] = findpeaks(Y, X, 'MinPeakHeight',Mean+5*Std, 'MinPeakDistance',1e3);
    relativeHeights = Heights - Mean;
    modeTable = array2table([relativeHeights, Heights, pksLocation], 'VariableNames', {'relativeHeights', 'Heights', 'pksLocation'});
end