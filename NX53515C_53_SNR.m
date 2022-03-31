%simulation parameters
global L w0x w0y lamda omegaM Gamma meff rho a b h bw P WtoV km k T kb hbar c SthFF Responsitivity PDGain
L = 55e-3;
w0x = 15e-6;
w0y = 15e-6;
lamda = 1064e-9;
omegaM = 122037*2*pi;
Gamma = 130.71;
rho = 3100;
a = 1.5e-3;
b = 3.5e-3;
h = 100e-9;
bw =8.515;
P = 1680e-6;
Responsitivity = 0.63;
PDGain = 1e4;
WtoV = Responsitivity*PDGain;
T = 300;
kb = 1.38064852e-23;
hbar = 1.054571817e-34;
c = 299792458;

meff = rho*a*b*h/4;
k = 2*pi/lamda;
km = 2*pi/a;
SthFF = 4*kb*T*Gamma*meff;

%simulation
wx = w(w0x,L);
ds = 0:wx/20:3*wx;
peakHeightC = zeros(1, length(ds));
noiFloorC = zeros(1, length(ds));
DN = -110;
samplepoints = [omegaM - 2*pi*5*bw:2*pi*bw:omegaM + 2*pi*5*bw];

for i = 1 : length(ds)
    d = ds(i);
    SPF = F(0,0,1,0,d);
    GeomNoise = GeomNoi(d, SPF);
    heights_C = arrayfun(@(x) S(x,pi/2, SPF, GeomNoise), samplepoints);
    power = sum(10.^(heights_C/10))/20;
    peakHeightC(i) = power;
    noiFloorC(i) = 10^(S(omegaM - 100*Gamma, pi/2, SPF, GeomNoise)/10)/20;
end
%peakHeightC = 10*log10(10.^(peakHeightC/10)+10^(DN/10));
%noiFloorC = 10*log10(10.^(noiFloorC/10)+10^(DN/10));

clear ans d f GeomNoise re SNR SPF wx y DN
ds = ds/w(w0x, L);

%datapart
Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\53-SNR\';

cPowerFlienameRegx = 'cPower_GS=*.bin';
cPowerFilelist = dir([Directory,cPowerFlienameRegx]);
[~,index] = sortrows({cPowerFilelist.date}.');
cPowerFilelist = cPowerFilelist(index); 
clear index;

powerFlienameRegx = 'power_GS=*.bin';
powerFilelist = dir([Directory,powerFlienameRegx]);
[~,index] = sortrows({powerFilelist.date}.');
powerFilelist = powerFilelist(index); 
clear index;

PSDoverP_FlienameRegx = 'PSDoverP_GS=*.bin';
PSDoverP_Filelist = dir([Directory,PSDoverP_FlienameRegx]);
[~,index] = sortrows({PSDoverP_Filelist.date}.');
PSDoverP_Filelist = PSDoverP_Filelist(index); 
clear index;

PSDoverPP_FlienameRegx = 'PSDoverPP_GS=*.bin';
PSDoverPP_Filelist = dir([Directory,PSDoverPP_FlienameRegx]);
[~,index] = sortrows({PSDoverPP_Filelist.date}.');
PSDoverPP_Filelist = PSDoverPP_Filelist(index); 
clear index;

PSDoverP_std_FlienameRegx = 'PSDoverP_std_GS=*.bin';
PSDoverP_std_Filelist = dir([Directory,PSDoverP_std_FlienameRegx]);
[~,index] = sortrows({PSDoverP_std_Filelist.date}.');
PSDoverP_std_Filelist = PSDoverP_std_Filelist(index); 
clear index;

PSDoverPP_std_FlienameRegx = 'PSDoverPP_std_GS=*.bin';
PSDoverPP_std_Filelist = dir([Directory,PSDoverPP_std_FlienameRegx]);
[~,index] = sortrows({PSDoverPP_std_Filelist.date}.');
PSDoverPP_std_Filelist = PSDoverPP_std_Filelist(index); 
clear index;

fileID = fopen([Directory, 'PSDfreq.bin']);
freq = fread(fileID, 'double');
fclose(fileID);

fileID = fopen([Directory, 'DN.bin']);
DN = fread(fileID, 'double');
fclose(fileID);

info = jsondecode(fileread([Directory,'info.json']));

peakF = 122037;
span = 500;
%peakF = 142980;
%span = 500;
noiFloorSampleRange = [113000, 116000];
%noiFloorSampleRange = [152000, 155000];

peakSearchRange = [peakF-span/2, peakF+span/2];
[~,PI1] = min(abs(freq-peakSearchRange(1)));
[~,PI2] = min(abs(freq-peakSearchRange(2)));

[~,NI1] = min(abs(freq-noiFloorSampleRange(1)));
[~,NI2] = min(abs(freq-noiFloorSampleRange(2)));

gapsizes = zeros(1, size(PSDoverP_Filelist,1));
steps = length(gapsizes);

for i =1:size(PSDoverP_Filelist, 1)
    temp = regexp(PSDoverP_Filelist(i).name, 'PSDoverP_GS=([0-9]+\.[0-9]*)\.bin', 'tokens');
    gapsizes(i) = str2num(temp{1}{1});
end

avg = 1000;
cPowers = zeros(1, steps);
for i =1:steps
    fileID = fopen([cPowerFilelist(i).folder,'\',cPowerFilelist(i).name]);
    data = fread(fileID, 'double');
    cPowers(i) = mean(data);
    fclose(fileID);
end

powers = zeros(1, steps);
n_powers = zeros(1, steps);
for i =1:steps
    fileID = fopen([powerFilelist(i).folder,'\',powerFilelist(i).name]);
    data = fread(fileID, 'double');
    fclose(fileID);
    
    fileID = fopen([cPowerFilelist(i).folder,'\',cPowerFilelist(i).name]);
    data2 = fread(fileID, 'double');
    fclose(fileID);
    
    powers(i) = mean(data);
    n_powers(i) = mean(data./data2)*mean(cPowers);
end


nois = zeros(1, length(gapsizes));
heights = zeros(1, length(gapsizes));
noi_stds = zeros(1, length(gapsizes));
height_stds = zeros(1, length(gapsizes));
SNRoverP = zeros(1, length(gapsizes));
for i = 1 : steps
    
    fileID = fopen([PSDoverP_Filelist(i).folder,'\',PSDoverP_Filelist(i).name]);
    data = fread(fileID, 'double');
    PSD = data *cPowers(i);
    fclose(fileID);
    PSD = PSD - DN;

    nois(i) = mean(PSD(NI1:NI2));
    
    fileID = fopen([PSDoverP_std_Filelist(i).folder,'\',PSDoverP_std_Filelist(i).name]);
    data = fread(fileID, 'double');
    std = data *cPowers(i);
    fclose(fileID);
    noi_stds(i) = sqrt(sum(std(NI1:NI2).^2)/(NI2-NI1+1)^2);
    
    fileID = fopen([PSDoverPP_Filelist(i).folder,'\',PSDoverPP_Filelist(i).name]);
    data = fread(fileID, 'double');
    PSD = data *cPowers(i)*cPowers(i);
    fclose(fileID);
    PSD = PSD - DN;
    PSD = PSD(PI1:PI2);
    [~,loc] = findpeaks(PSD, 'MinPeakDistance', PI2-PI1-1);
    heights(i) = sum(PSD(loc-5:loc+5));
    
    fileID = fopen([PSDoverPP_std_Filelist(i).folder,'\',PSDoverPP_std_Filelist(i).name]);
    data = fread(fileID, 'double');
    std = data *cPowers(i)*cPowers(i);
    std = std(PI1:PI2);
    fclose(fileID);
    height_stds(i) = sqrt(sum(std(loc-5:loc+5).^2));
end
SNR = heights./nois;
SNR_std = sqrt((1./nois).^2.*height_stds.^2 + (heights./nois.^2).^2.*noi_stds.^2);
gapsizes = max(gapsizes)-gapsizes;

% fit
fitMirror = 'A';
if fitMirror == 'A'
    Axfit = erfFit(info.MirrorA.positions.value, n_powers, info.MirrorA.angle.value, info.Mirror_tilt_angle.value, 7.3*1.5*1.64, fitMirror);
    GSfit = 2*(Axfit.x0*1e3-info.MirrorA.positions.value)*sin(info.MirrorA.angle.value)*cos(info.Mirror_tilt_angle.value);
    %GSfit = 2*(8.67-info.MirrorA.positions.value)*sin(info.MirrorA.angle.value)*cos(info.Mirror_tilt_angle.value);
    GSfit_overW = GSfit/Axfit.wx*1e-3;
else
    Bxfit = erfFit(info.MirrorB.positions.value, n_powers, info.MirrorB.angle.value, info.Mirror_tilt_angle.value, 7.3*1.5*1.64, fitMirror);
    GSfit = 2*(Bxfit.x0*1e3-info.MirrorB.positions.value)*sin(info.MirrorB.angle.value)*cos(info.Mirror_tilt_angle.value);
    GSfit_overW = GSfit/Bxfit.wx*1e-3;
end

%plot
[~,cut] = min(abs(ds-GSfit_overW(end)));
%cut = cut+1;
noi_offset = -6.3;
height_offset = -3;

subplot(2,3,1);
height_logstds = 10*height_stds./(log(10)*heights);
errorbar(GSfit_overW, 10*log10(20*heights), height_logstds,'.');
hold on;
plot(ds(1:cut), 10*log10(20*peakHeightC(1:cut)));
title('Subplot 1: Peak area');
xlabel('Gapsize/beamWaist')
ylabel('Power(mdB)')

subplot(2,3,2);
noi_logstds = 10*noi_stds./(log(10)*nois);
errorbar(GSfit_overW, 10*log10(20*nois), noi_logstds, '.');
hold on;
plot(ds(1:cut), 10*log10(20*noiFloorC(1:cut)));
title('Subplot 2: Shot noise');
xlabel('Gapsize/beamWaist')
ylabel('Power(mdB)')

subplot(2,3,3);
errorbar(GSfit_overW, 10*log10(20*SNR), sqrt(noi_logstds.^2+height_logstds.^2),'.');
hold on;
plot(ds(1:cut), 10*log10(20*peakHeightC(1:cut)./noiFloorC(1:cut)));
title('Subplot 3: SNR');
xlabel('Gapsize/beamWaist')
ylabel('Power(mdB)')

subplot(2,3,4);
height_logstds = 10*height_stds./(log(10)*heights);
errorbar(GSfit_overW, 10*log10(20*heights)+height_offset, height_logstds,'.');
hold on;
plot(ds(1:cut), 10*log10(20*peakHeightC(1:cut)));
title(sprintf('Subplot 4: Peak area, offseted by %0.2f dB', height_offset));
xlabel('Gapsize/beamWaist')
ylabel('Power(mdB)')

subplot(2,3,5);
noi_logstds = 10*noi_stds./(log(10)*nois);
errorbar(GSfit_overW, 10*log10(20*nois)+noi_offset, noi_logstds, '.');
hold on;
plot(ds(1:cut), 10*log10(20*noiFloorC(1:cut)));
title(sprintf('Subplot 5: Shot Noise, offseted by %0.2f dB', noi_offset));
xlabel('Gapsize/beamWaist')
ylabel('Power(mdB)')

subplot(2,3,6);
errorbar(GSfit_overW, 10*log10(20*SNR)+height_offset-noi_offset, sqrt(noi_logstds.^2+height_logstds.^2),'.');
hold on;
plot(ds(1:cut), 10*log10(20*peakHeightC(1:cut)./noiFloorC(1:cut)));
title(sprintf('Subplot 6: SNR, offseted by %0.2f dB', height_offset-noi_offset));
xlabel('Gapsize/beamWaist')
ylabel('Power(mdB)')

set(gcf, 'Position',  [100, 100, 1350, 900])

%functions
function re = SNRc(omega, theta, d)
    global Gamma
    re = S(oemga, theta, d) - S(oemga-100*Gamma, theta, d);
end

function re = chi(omega)
    global meff omegaM Gamma
    re = 1/meff*((omegaM^2 - omega.^2)+1i*Gamma*omega).^(-1);
end

function re = GeomNoi(d, denominator)
    re = 0;
    for i = 1:5
        re = re+(abs(F(0,0,2*i+1,0,d)/denominator)).^2;
    end
end

function re = F(m,n,l,k,d)
    global w0x w0y L
    wx = w(w0x, L);
    wy = w(w0y, L);
    scale = 3;
    func = @(x, y) u(m,l,x,x).*Fweight(x,d);
    re = integral(func, -scale*wx, scale*wx,'RelTol',1e-3);
end

function re = Fweight(x,d)
    re = heaviside(x - d/2) - heaviside(-(x + d/2));
end

function re = u(m,n,x,y)
    global w0x w0y L
    wx = w(w0x, L);
    wy = w(w0y, L);
    re = (2/pi)^(1/2)/sqrt(2^(m+n)*factorial(m)*factorial(n)*wx*wy)*exp(-x.^2/wx^2-y.^2/wy^2).*hermiteH(m,sqrt(2)*x/wx).*hermiteH(n,sqrt(2)*y/wy);
end

function re = w(w0, z)
    global lamda
    re = w0*sqrt(1+(z*lamda/(pi*w0^2))^2);
end

function re = S(omega, theta, SPF, GeomNoise)
    global P bw k km w0x SthFF WtoV hbar c Responsitivity
    re = 10*log10(20) + ...
    10*log10(2*bw*((2*P*SPF*k*km*w0x*abs(chi(omega)).*sin(theta)).^2*SthFF/2 +...
    P*SPF^2*hbar*c*k*(abs(2*P/(c*k)*(k*km*w0x)^2*chi(omega).*sin(theta) - cos(theta)).^2/Responsitivity + sin(theta).^2/Responsitivity + GeomNoise/Responsitivity))*(WtoV)^2); %in unit of dB/bin
end

function re = Sth(f, theta, SPF, GeomNoise)
    global P bw k km w0x SthFF WtoV hbar c
    re = 10*log10(20*bw*4*P^2*SPF^2*(k*km*w0x)^2*SthFF*(abs(chi(2*pi*f))).^2*(WtoV)^2);
end

omega = c*2*pi/lamda;
P/(omega*meff*omegaM*Gamma)*k^2*km^2*w0x^2