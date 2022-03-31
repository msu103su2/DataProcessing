%experiment data
Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\manuscript_data\';
fileID = fopen([Directory, 'psd_avg20_50s.bin']);
PSD = fread(fileID, 'double');
fclose(fileID);

%abotu 90mW on device, 54mW thransmitted, 33mW reflected

fileID = fopen([Directory, 'freqs_50s.bin']);
freqs = fread(fileID, 'double');
fclose(fileID);
    
%device properties
rho = 3100;
b = 1.5e-3;
a = 3.5e-3;
h = 100e-9;

%mode properties

fm = 123774;
gamma = 24.4074;
modeIndex = 4;

fm = 223339;
gamma = 6.5481;
modeIndex = 8;

fm = 275324;
gamma = 2.2662;
modeIndex = 10;

%fm = 326472.8;
%gamma = 1.144;
%modeIndex = 10;

fm = 371696.906;
gamma = 0.7295;
modeIndex = 12;
%laser properties
w0 = 165e-6;
lenf = 125e-3;
lambda = 1064e-9;
Pdc = 27.2e-3;

%chirp properties
chirpSpanf = 100;
chirp_fc = 218533.9;
DeltaX0 = 7.8e-6;

%photodetector properties
WtoA = 0.75;
pdGain = 2/3*100*1*1e3;
gapsize =18e-6;

%DAQ properties
fs = 1e7;
resolution = 0.1;%Hz

%circuit attenuation:
attenuation = 0.88*0.5623; %in voltage

%set experiment
sp = SplitPair(gapsize, 0);
device = Membrane(a, b, h, rho);
device.SetMode(fm, gamma, modeIndex);
laser = Laser(lambda);
laser.SetPdc(Pdc);
laser.deltaf = chirpSpanf;
PD = PhotoDetector(WtoA, pdGain);
f = device.fm-50000*device.gamma:0.01:device.fm+50000*device.gamma;
OLE = OLExperiment(sp, device, laser, w0, lenf, PD);
OLE.SetInspectFreqs(f);
OLE.CalF(0, gapsize);

%calculation
[total, signal, noise, thermal, Qnoi, ref] = OLE.GetPSDs_AOM(1.55, DeltaX0);

%apply circuit attenuation
total = total + 2*10*log10(attenuation);
signal = signal + 2*10*log10(attenuation);
noise = noise + 2*10*log10(attenuation);
thermal = thermal + 2*10*log10(attenuation);
Qnoi = Qnoi + 2*10*log10(attenuation);

factor = 1/(2*pi);
check = 10.^(thermal/10)*factor + 10.^(Qnoi/10);
check = 10*log10(check);

%calculate C=1
hbar = 1.0545718e-34;
kb = 1.380649e-23;
T = 273.15+24.37;
C_1 = 10*log10(kb*T./(hbar*2*pi.*f)) + mean(Qnoi);

%convert back to sxx
elec = 2*10*log10(attenuation);
G = 4*Pdc^2*OLE.F1^2*(w0*OLE.device.kmx()*OLE.laser.k)^2;
G = 10*log10(G);

%plot
figure()
fit = plot(f, check, 'b', 'DisplayName','calculation');
hold on
data = plot(freqs, PSD, '.', 'DisplayName','data');
C_1_plot = plot(f, C_1, 'g-', 'DisplayName','C=1');
uistack(fit,'top')
xlabel('f(Hz)')
ylabel('power(dBm/Hz)')
legend()

%plotSxx
figure()
fit = semilogy(f, 10.^((check - elec-G-10*log10(20))/10), 'b', 'DisplayName','calculation');
hold on
data = semilogy(freqs, 10.^((PSD - elec-G-10*log10(20))/10), '.', 'DisplayName','data');
C_1_plot = semilogy(f, 10.^((C_1 - elec-G-10*log10(20))/10), 'g-', 'DisplayName','C=1');
uistack(fit,'top')
xlabel('f(Hz)')
ylabel('Sxx(m^2/Hz)')
legend()


%application
phis = 1.561 : 0.001/5:1.581;
Z = zeros(length(f), length(phis));
i = 1;
for phi = phis
    [total, signal, noise, thermal, ref] = OLE.GetPSDs_AOM(phi, DeltaX0);
    Z(:,i) = signal - ref;
    i = i+1;
end
Z = Z.';

idx = [1:49,51:101];

Y = pi-phis;
X = f - fm;
figure
[C,h] = contourf(X, Y(idx), Z(idx, :), 50);
set(h,'LineColor','none');
colormap(jet)