function FindPeaks(X,Y)
    [pks,f] = findpeaks(Y, X, 'MinPeakHeight',-95, 'MinPeakDistance',1e3);
    plot(1,f,'*');
end