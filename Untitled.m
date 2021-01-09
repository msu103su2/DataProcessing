for i = 1:100
    data = importdata([filelist(i).folder,'\',filelist(i).name]);
    data = data.data;
    checkrange = [121000, 123000];
    [~,i1] = min(abs(f-checkrange(1)));
    [~,i2] = min(abs(f-checkrange(2)));
    VA = data(:,2)*1;
    VB = data(:,3);
    VC = data(:,4);
    [psd, f] = periodogram(VA-VB, rectwin(length(VA)), length(VA), 1/(data(2,1)-data(1,1)));
    [psdC, f] = periodogram(VC, rectwin(length(VC)), length(VC), 1/(data(2,1)-data(1,1)));
    [~,i1] = min(abs(f-checkrange(1)));
    [~,i2] = min(abs(f-checkrange(2)));
    plot(f(i1:i2),10*log10(psd(i1:i2)*20), f(i1:i2), 10*log10(psdC(i1:i2)*20));
    waitforbuttonpress
end