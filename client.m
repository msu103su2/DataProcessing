folder = 'Z:/data/optical lever project/Die 000_21/03_scanPSDs/';
WaferSN = '000';
DieNumber = 21;
fileIndex = [1:14,16,17];
FrequencyRange = [1e5,1e6];
filenameRegxs = [];
global samplesPerFile;
samplesPerFile = 801;

deviceid = [1,2,3,4,5,6,0,7,8,11,12,0,13,14,15,16,17];
warning('off')
ModesofAllDevice = {};
figure('Name','WaterFall');
hold on;
for i = fileIndex
    data = concatenate(folder, sprintf('%s_%i_%02i_*.csv',WaferSN, DieNumber, i),FrequencyRange);
    index = data.FrequenciesHz < FrequencyRange(2) & data.FrequenciesHz > FrequencyRange(1);
    if (deviceid(i)~=0)
        plot(data.PSDdB2Hz(index)+deviceid(i)*100,data.FrequenciesHz(index)/1e6);
    end
    modeTable = table('size', [1 3], 'VariableNames', {'relativeHeights', 'Heights', 'pksLocation'}, 'VariableTypes',{'double','double','double'});
    
    for j = 1 : samplesPerFile : size(data.FrequenciesHz,1)
        modeTable = [modeTable; FindPeaks(table2array(data(j : min(j+samplesPerFile-1, size(data.FrequenciesHz,1)), 'FrequenciesHz')), table2array(data(j : min(j+samplesPerFile-1, size(data.FrequenciesHz,1)), 'PSDdB2Hz')))];
    end
    
    modeTable(1,:) = [];
    modeTable = sortrows(modeTable, 3);

    ModesofAllDevice{end+1} = modeTable;
end
hold off;


mark = [];
for i = 1 : size(ModesofAllDevice,2)
    if(isempty(ModesofAllDevice{i}))
        mark = [mark, i];
    end
end
ModesofAllDevice(mark) = [];

figure('Name','Measured Modes');
hold on;
for i = 1 : size(ModesofAllDevice,2)
    indexs = putTogether(ModesofAllDevice{i}.pksLocation, 1e4);
    if (deviceid(i)~=0)
        plot3(zeros(size(indexs,2),1)+deviceid(i),ModesofAllDevice{i}.pksLocation(indexs),ModesofAllDevice{i}.relativeHeights(indexs),'*');
    end
end
hold off;


%{
sz = 0;
for i = 1 : size(ModesofAllDevice,2)
    sz = sz + size(ModesofAllDevice{i},1);
end

AllModesTable = table('size', [sz 4], 'VariableNames', {'relativeHeights', 'Heights', 'pksLocation', 'DeviceNumber'}, 'VariableTypes',{'double','double','double', 'int16'});
index = 1;
for i = 1 : size(ModesofAllDevice,2)
    AllModesTable(index: index+size(ModesofAllDevice{i},1)-1,{'relativeHeights', 'Heights', 'pksLocation'}) = ModesofAllDevice{i}(:,{'relativeHeights', 'Heights', 'pksLocation'});
    AllModesTable.DeviceNumber(index:index+size(ModesofAllDevice{i},1)-1) = zeros(size(ModesofAllDevice{i},1),1)+i;
    index = index + size(ModesofAllDevice{i},1);
end
AllModesTable = sortrows(AllModesTable, 'pksLocation');
i = 1;
ModeFs = [];
counts = [];
while i < size(AllModesTable,1)
    j = i + 1;
    count = 1;
    while j < size(AllModesTable,1) + 1
        if (issame(AllModesTable.pksLocation(i), AllModesTable.pksLocation(j)))
            j = j + 1;
            count = count + 1;
        else
            i = j;
            break;
        end
    end
    ModeFs = [ModeFs; AllModesTable.pksLocation(i)];
    counts = [counts; count];
end
commonModes = table('size',[size(ModeFs,1) 2],'VariableNames', {'ModeF', 'count'}, 'VariableTypes',{'double', 'int16'});
commonModes.ModeF = ModeFs;
commonModes.count = counts;
figure('Name','Common modes')
plot(commonModes.ModeF,commonModes.count, '*');
%}
warning('on')

function modeTable = FindPeaks(X,Y)
    Mean = mean(Y);
    Std = std(Y);
    [Heights,pksLocation] = findpeaks(Y, X, 'MinPeakHeight',Mean+2*Std, 'MinPeakDistance',1e3);
    relativeHeights = Heights - Mean;
    modeTable = array2table([relativeHeights, Heights, pksLocation], 'VariableNames', {'relativeHeights', 'Heights', 'pksLocation'});
end

function dataTable = concatenate(folder, flienameRegx, FrequencyRange)
    global samplesPerFile;
    filelist = dir([folder,flienameRegx]);
    dataTable = table('Size',[samplesPerFile*size(filelist(1).name,1) 2],'VariableNames', {'FrequenciesHz', 'PSDdB2Hz'},'VariableTypes', {'double', 'double'});
    for i = 1 : size(filelist,1)
        newtable = importfile([folder,filelist(i).name]);
        dataTable((i-1)*samplesPerFile + 1 : i * samplesPerFile ,{'FrequenciesHz', 'PSDdB2Hz'}) = newtable(:,{'Frequencies', 'PSDdB2Hz'});
    end
    dataTable = dataTable(dataTable.FrequenciesHz > FrequencyRange(1) & dataTable.FrequenciesHz < FrequencyRange(2),:);
end

function result = issame(fx,fy)
    if abs(fx-fy) < 1e3
        result = true;
    else
        result = false;
    end
end

function indexs = putTogether(pksLocation, frange)
    i = 1;
    indexs = [];
    while i < size(pksLocation, 1)-1
        j = i + 1;
        while j < size(pksLocation, 1) + 1
            if pksLocation(j) - pksLocation(i) < frange && j < size(pksLocation, 1)
                j = j + 1;
            else
                avg = mean(pksLocation(i:j-1));
                [temp,index] = min(abs(pksLocation(i:j-1) - avg));
                indexs = [indexs, index + i - 1];
                i = j;
                if (i < size(pksLocation, 1)-1)
                    break;
                else
                    return;
                end
            end
        end
    end
end