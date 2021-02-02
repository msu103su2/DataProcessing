function pbest = PDdualAC_balance2(data)
    VC = data(:,3);
    rate = 1e8/6.4;
    inspectF = 122000;
    span = 20000;

    [psdC,f] = periodogram(VC, rectwin(length(VC)), length(VC), rate);
    [~,i1] = min(abs(f-inspectF+span));
    [~,i2] = min(abs(f-inspectF-span));

    p = 10*log10(psdC(i1:i2)*20);
    [~, loc] = findpeaks(p, 'MinPeakProminence',50);
    
    VA = data(:,1);
    VB = data(:,2);
    
    fftA = fft(VA);
    fftB = fft(VB);
    N = length(VA);
    
    if (rem(N,2)==0)
        fftA = fftA(1:N/2+1);
        fftB = fftB(1:N/2+1);
    else
        fftA = fftA(1:(N+1)/2);
        fftB = fftB(1:(N+1)/2);
    end
    
    data = [fftA(loc+i1-1),fftB(loc+i1-1)];

    funcO = @(x) func(data, x);
    pbest = fminbnd(funcO, 0.9, 1.1);
end

function meritValue = func(data, p)
    VA = data(:,1)*p;
    VB = data(:,2)*(2-p);
    meritValue = sum((abs(VA-VB)).^2);
end