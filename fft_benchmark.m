file = 'Z:\data\optical lever project\NORCADA_NX53515C\20-SNR\TS_GS=12.300440899753948_count=999.bin';
fileID = fopen(file);
data = (fread(fileID, [3,inf], 'double'))';
VA = [data(:,1);0];
VB = [data(:,2);0];

avg = 1;
rate = 1e8/6.4;

PDdualAC_balance2(data)

%{
[psdref,fref] = periodogram(VA-VB, rectwin(length(VA)), length(VA), rate);
Sk = fft(VA-VB);

N = length(VA);
T = N/rate;
if (rem(N,2)==0)
    f = [0:N/2]/T;
    Sk = Sk(1:N/2+1);
else
    f = [0:(N-1)/2]/T;
    Sk = Sk(1:(N+1)/2);
end

psd = 2*abs(Sk).^2/(T*rate*rate);
psd(1) = psd(1)/2;

if (rem(N,2)==0)
    psd(end) = psd(end)/2;
end

psd1 = 10*log10(20*psd);
psd2 = 10*log10(20*psdref);
plot(1:10,psd1(1:10),1:10, psd2(1:10))
%}
