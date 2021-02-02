function pbest = PDdualAC_balance(data)
    VC = data(:,3);
    rate = 1e8/6.4;
    inspectF = 122000;
    span = 20000;

    [psdC,f] = periodogram(VC, rectwin(length(VC)), length(VC), rate);
    [~,i1] = min(abs(f-inspectF+span));
    [~,i2] = min(abs(f-inspectF-span));

    p = 10*log10(psdC(i1:i2)*20);
    f = f(i1:i2);
    [pks, loc] = findpeaks(p, 'MinPeakProminence',50);

    funcO = @(x) func(data, x, f(loc), rate);
    pbest = fminbnd(funcO, 0.9, 1.1);
end

function meritValue = func(data, p, f, rate)
    VA = data(:,1)*p;
    VB = data(:,2)*(2-p);
    [peaks, ~] = periodogram(VA-VB, rectwin(length(VA)), f, rate);
    meritValue = 10*log10(sum(peaks)*20);
end