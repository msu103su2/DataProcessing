data = importdata([filelist(1).folder,'\',filelist(1).name]);
data = data.data;
checkrange = [12000, 500000];
[~,i1] = min(abs(f-checkrange(1)));
[~,i2] = min(abs(f-checkrange(2)));
for p = 0.98:0.005:1.03
    VA = data(:,2);
    VB = data(:,3);
    %VC = data(:,4);
    [psd, f] = periodogram(VA-VB, rectwin(length(VA)), length(VA), 1/(data(2,1)-data(1,1)));
    [psdC, f] = periodogram(VC, rectwin(length(VC)), length(VC), 1/(data(2,1)-data(1,1)));
    [~,i1] = min(abs(f-checkrange(1)));
    [~,i2] = min(abs(f-checkrange(2)));
    plot(f(i1:i2),10*log10(psd(i1:i2)*20));
    waitforbuttonpress;
end