folder = 'Z:/data/optical lever project/NORCADA_NX53515C/02-balancing';
global samplesPerFile;
samplesPerFile = 801;
OPmax = 177.1829;
transmitted = [127.28, 125.77, 124.51,123.92, 124.03,124.18,124.92,126.32];
OP = OPmax - transmitted;

%{
for  j = 22:2:36
    flienameRegx = ['000_13_test_' , num2str(j),'*.csv']
    filelist = dir([folder,flienameRegx]);
    dataTable = table('Size',[samplesPerFile*size(filelist(1).name,1) 2],'VariableNames', {'FrequenciesHz', 'PSDdB2Hz'},'VariableTypes', {'double', 'double'});
    for i = 1 : size(filelist,1)
        newtable = importfile1([folder,filelist(i).name]);
        dataTable((i-1)*samplesPerFile + 1 : i * samplesPerFile ,{'FrequenciesHz', 'PSDdB2Hz'}) = newtable(:,{'FrequenciesHz', 'PSDdB2Hz'});
    end
    FindPeaks(dataTable.FrequenciesHz, dataTable.PSDdB2Hz)
    figure('Name',num2str(j));
    plot(dataTable.FrequenciesHz, dataTable.PSDdB2Hz);
end
%}
warning('off')
pksLocation = [];
heights = [];
position = [];
positioni = 1;
for  j = 36:-2:22
    flienameRegx = '*.csv'
    filelist = dir([folder,flienameRegx]);
    for i = 1 : size(filelist,1)
        newTable = importfile1([folder,filelist(i).name]);
        dataTable = FindPeaks(newTable.FrequenciesHz, newTable.PSDdB2Hz);
        pksLocation = [pksLocation; dataTable.pksLocation];
        heights = [heights; dataTable.relativeHeights - 10*log10(OP(positioni))];
        position = [position; zeros(size(dataTable,1),1)+j];
    end
    positioni = positioni+1;
    %figure('Name',num2str(j));
    %plot(dataTable.FrequenciesHz, dataTable.PSDdB2Hz);
end
[fMean,fStd] = group(300, position, pksLocation, heights);
fileID = fopen('Z:/data/optical lever project/pull/modes.txt','w');
nbytes = fprintf(fileID,'%d %d\n',[fMean,fStd]');
fclose(fileID);
warning('on')
%plot3(position, pksLocation, heights,'*');

function modeTable = FindPeaks(X,Y)
    Mean = mean(Y);
    Std = std(Y);
    [Heights,pksLocation] = findpeaks(Y, X, 'MinPeakHeight',Mean+3*Std, 'MinPeakDistance',1e3);
    relativeHeights = Heights - Mean;
    modeTable = array2table([relativeHeights, Heights, pksLocation], 'VariableNames', {'relativeHeights', 'Heights', 'pksLocation'});
end

function [fMean,fStd] = group(width, position, pksLocation, heights)
    [pksLocation, order] = sort(pksLocation);
    heights = heights(order);
    position = position(order);
    i = 1;
    fMean = [];
    fStd = [];
    hold on;
    while i < size(pksLocation, 1)
        for j = (i+1):size(pksLocation, 1)
            if pksLocation(j) - pksLocation(i) > width
                [x, order] = sort(position(i : j-1));
                z = heights(i : j-1);
                y = pksLocation(i : j-1);
                z = z(order);
                y = y(order);
                fMean = [fMean; mean(y)];
                fStd = [fStd; std(y)];
                plot3(x, y, z);
                i = j;
                break;
            end
            
            if j == size(pksLocation, 1)
                    i = j;
            end
        end
    end
end