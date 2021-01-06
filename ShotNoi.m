Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\03-ModeCharacterizaiton\';
flienameRegx = 'BP*R.csv'
filelist = dir([Directory,flienameRegx]);

noiFloorSampleRange = [112000, 116000];
noiFloor = [];
opticalPower = [];

data = importdata([Directory,'BPDNR.csv']);
data = data.data;
data(:,1) = data(:,1)*1e6;
f = data(:,1);
p = data(:,2);
[~,Indexs(1)] = min(abs(f-noiFloorSampleRange(1)));
[~,Indexs(2)] = min(abs(f-noiFloorSampleRange(2)));
darkNoi = 10.^(p(Indexs(1):Indexs(2))/10)/20;

for i=1:size(filelist,1)
    if(strcmp(filelist(i).name,'BPDNR.csv'))
    else
        data = importdata([filelist(i).folder,'\',filelist(i).name]);
        data = data.data;
        data(:,1) = data(:,1)*1e6;

        f = data(:,1);
        p = data(:,2);

        [~,Indexs(1)] = min(abs(f-noiFloorSampleRange(1)));
        [~,Indexs(2)] = min(abs(f-noiFloorSampleRange(2)));

        noiFloor = [noiFloor, mean(10.^(p(Indexs(1):Indexs(2))/10)/20 - darkNoi)];

        match = regexp(filelist(i).name,'BP([0-9]+.*[0-9]*)R','tokens');
        opticalPower = [opticalPower, str2double(match{1}{1})];
    end
end

x = opticalPower;
y = noiFloor;
p = polyfit(x, y, 1);
yfit = polyval(p, x);
SSresid = sum((y - yfit).^2);
SStotal = (length(y)-1)*var(y);
rsq = 1 - SSresid/SStotal;

figure(1);
hold on;
plot(opticalPower, noiFloor,'*');
fplot(@(X) p(1)*X+p(2));
axis([0 2500 0 1e-11]);
legend('Data','LinearFit');
xlabel('LaserPower(uW)');
ylabel('BalancedNoisePower(W)');
dim = [.2 .5 .3 .3];
annotation('textbox',dim,'String','R2 = 0.995','FitBoxToText','on');