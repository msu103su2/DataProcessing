Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\167-characterize\';

PSDFlienameRegx = 'PSD_*.bin';
PSDFilelist = dir([Directory,PSDFlienameRegx]);
[~,index] = sortrows({PSDFilelist.date}.');
PSDFilelist = PSDFilelist(index); 
clear index;

DrivePSDFlienameRegx = 'DrivePSD_*.bin';
Filelist = dir([Directory,DrivePSDFlienameRegx]);
[~,index] = sortrows({Filelist.date}.');
Filelist = Filelist(index); 
clear index;

fileID = fopen([Directory, 'freq.bin']);
freq = fread(fileID, 'double');
fclose(fileID);

p = zeros(1,length(PSDFilelist));

for i = 1:length(PSDFilelist)
    fileID = fopen([PSDFilelist(i).folder,'\',PSDFilelist(i).name]);
    PSD = fread(fileID, 'double');
    fclose(fileID);
    [p(i), I] = max(PSD);
    plot(freq, PSD);
    waitforbuttonpress;
    fileID = fopen([Filelist(i).folder,'\',Filelist(i).name]);
    DrivePSD = fread(fileID, 'double');
    fclose(fileID);
    p(i) = p(i)-DrivePSD(I);
end
plot(p)

SPoffset_start = -60e-6;
SPoffset_step = 6e-6;
SPoffset_end = 60e-6;

laseroffset_start = -60e-6;
laseroffset_step = 6e-6;
laseroffset_end = 60e-6;

result = zeros(length(SPoffset_start:SPoffset_step:SPoffset_end), length(laseroffset_start:laseroffset_step:laseroffset_end));
index = 1;
for SPoffset = SPoffset_start:SPoffset_step:SPoffset_end
    SP.offset = SPoffset;
    F0 = SP.F_onlyx(HGB0, 0);
    F1 = SP.F_onlyx(HGB1, 0);
    laseroffset = laseroffset_start:laseroffset_step:laseroffset_end;
    result(index, :) = F0 + laseroffset/60e-6*F1;
    index = index+1;
end
%%%
rho = 3100;
b = 1.5e-3;
a = 3.5e-3;
h = 100e-9;
w0 = 120e-6;
lenf = 125e-3;
lambda = 1064e-9;

fm = 218600;
gamma = 3;
modeIndex = 8;

Pdc = 4.3e-3;
chirpAC = 0.087e-3;
chirpSpanf = 100;
sidebandAC = 0.38e-3;

pdGain = 10000;

gapsize = 45e-6;
sp = SplitPair(gapsize, 0);
device = Membrane(a, b, h, rho);
device.SetMode(fm, gamma, modeIndex);
laser = Laser(lambda);
laser.SetPdc(Pdc);
laser.SetAC(chirpAC, chirpSpanf);
laser.Pac_sideband = sidebandAC;
PD = PhotoDetector(pdGain);
angle = HGBeam.guoyAngle(w0, lambda, flen);

f = device.fm-500*device.gamma:0.25:device.fm+500*device.gamma;
OLE.SetInspectFreqs(f);

DeviceMove = -20e-6;
DSMove = -0.5e-6;
mirrorGapSize = 100e-6;
laseroff = -DeviceMove;
DSoff = DSMove-DeviceMove;

OLE.CalF(DSoff, mirrorGapSize);
[total, signal, noise, signal_chi, thermal, sideband] = OLE.GetPSDs(laseroff, angle);
plot(x,y,f, noise);
plot(f, signal, f, noise);
%OLE.CoefA(laseroff, angle)
mean(OLE.DC_CM(laseroff, angle))
plot(f- OLE.device.fm, total, '*');

%plot(f, noise, ThCalifreq, ThCali);
plot(f - OLE.device.fm, twof, f- OLE.device.fm, signal, f- OLE.device.fm, noise);
mean(OLE.DC_CM(laseroff, angle))
%%%

spdxs = 0:0.01e-5:2e-5;
distance = zeros(1,length(spdxs));
[total, signal, noise] = OLE.GetPSDs(2e-5, angle);
plot(f, signal, f, noise);
for i = 1 : length(spdxs)
    OLE.CalF(spdxs(i), 100e-6);
    OLE.F0 = 0;
    [total, signal, noise] = OLE.GetPSDs(1e-5, angle);
    plot(f, signal);
    waitforbuttonpress;
    [~,minI] = min(signal);
    [~,maxI] = max(signal);
    distance(i) = f(maxI) - f(minI);
end
plot(spdxs*1e6, distance);
xlabel('SPoffset(um)');
ylabel('PPdf(Hz)');




Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\114-LineScan\';
PSDFlienameRegx = 'signal_*.bin';
PSDFilelist = dir([Directory,PSDFlienameRegx]);
[~,index] = sortrows({PSDFilelist.date}.');
PSDFilelist = PSDFilelist(index); 
clear index;

Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\138-LineScan\';
driveFlienameRegx = 'drive_*.bin';
driveFilelist = dir([Directory,driveFlienameRegx]);
[~,index] = sortrows({driveFilelist.date}.');
driveFilelist = driveFilelist(index); 
clear index;

fileID = fopen([Directory, 'freq.bin']);
freq = fread(fileID, 'double');
fclose(fileID);

f = device.modes(29.7);
inspectfs = f(39, 1);
span = 2e3;
flos = inspectfs - span/2;
fhis = inspectfs + span/2;
ilos = zeros(length(flos));
ihis = zeros(length(fhis));

for i = 1:length(inspectfs)
    [~,ilos(i)] = min(abs(freq - flos(i)));
    [~,ihis(i)] = min(abs(freq - fhis(i)));
end

data = zeros(length(inspectfs), length(PSDFilelist));
sumspan = 15;
for i = 1:length(PSDFilelist)
    fileID = fopen([Directory, PSDFilelist(i).name]);
    PSD = fread(fileID, 'double');
    fclose(fileID);

    fileID = fopen([Directory, driveFilelist(i).name]);
    drive = fread(fileID, 'double');
    fclose(fileID);
    
   
    for j = 1 : length(inspectfs)
        PSDj = PSD(ilos(j):ihis(j));
        drivej = drive(ilos(j):ihis(j));
        PSDj = 10*log10(PSDj);
        drivej = 10*log10(drivej);
        fj = freq(ilos(j):ihis(j));
        [~,maxI] = max(PSDj);
        %plot(freq(ilos(j):ihis(j)),10*log10(20*PSD(ilos(j):ihis(j))),freq(ilos(j):ihis(j)), 10*log10(20*drive(ilos(j):ihis(j))) );
        %waitforbuttonpress;
        temp = fit_logLorentzian(fj(maxI-300:maxI+300), PSDj(maxI-300:maxI+300) - drivej(maxI-300:maxI+300), fj(maxI));
        data(j,i) = temp.offset;
    end
end

fileID = fopen([Directory, 'DriveDFT_10120.bin']);
data = fread(fileID, 'double');
fclose(fileID);

data = complex(data(1:2:end), data(2:2:end));
temp = zeros(10, 1600);
for i = 1:10
    temp(i,:) = data(1600*(i-1)+1:1600*i);
end

data = zeros(10, 1600);
for i = 1:10
    data(i,:) = A(i,:)./B(i,:);
end
    
temp = zeros(10, 241);
for i = 1:10
    [~, I] = max(abs(A(i,:)));
    temp(i,:) = data(i, I-120:I+120);
end

Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\142-Dtest\';
fileID = fopen([Directory, 'DC_on_step.bin']);
DC_A = fread(fileID, 'double');
fclose(fileID);
fileID = fopen([Directory, 'powers.bin']);
powers = fread(fileID, 'double');
fclose(fileID);
fileID = fopen([Directory, 'spoffs.bin']);
spoffs = fread(fileID, 'double');
fclose(fileID);
fileID = fopen([Directory, 'spdxs.bin']);
spdxs = fread(fileID, 'double');
fclose(fileID);
fileID = fopen([Directory, 'Axs.bin']);
Axs = fread(fileID, 'double');
fclose(fileID);
fileID = fopen([Directory, 'Bxs.bin']);
Bxs = fread(fileID, 'double');
fclose(fileID);

Cali_A0 = 8.615;
Cali_Aw0 = 1.180;
Cali_B0 = 9.832;
Cali_Bw0 = 1.113;

Adxs = (Cali_A0 - Axs)/Cali_Aw0*HGB0.wx(OLE.lenf);
Bdxs = (Cali_B0 - Bxs)/Cali_Bw0*HGB0.wx(OLE.lenf);
spdxs_cali = Adxs + Bdxs;
spoffs_cali = (-Adxs + Bdxs)/2;
[~,I] = min(DC_A.^2);
spoffs_cali = spoffs_cali - spoffs_cali(I); 

chirp_span  = 20;
N = length(powers);
powers = powers*chirp_span;
[DC, AC] = Match_DS(OLE, spdxs_cali, spoffs_cali);
plot(spoffs_cali, -DC, spoffs_cali, DC_A, 'o')
plot(spoffs_cali, AC, spoffs_cali, powers, 'o')

Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\136-test\';
fileID = fopen([Directory, 'DC_on_step.bin']);
DC_A = fread(fileID, 'double');
fclose(fileID);
fileID = fopen([Directory, 'powers.bin']);
powers = fread(fileID, 'double');
fclose(fileID);
fileID = fopen([Directory, 'ref.bin']);
ref = fread(fileID, 'double');
fclose(fileID);
fileID = fopen([Directory, 'Axs.bin']);
Axs = fread(fileID, 'double');
fclose(fileID);
fileID = fopen([Directory, 'Bxs.bin']);
Bxs = fread(fileID, 'double');
fclose(fileID);
fileID = fopen([Directory, 'freq.bin']);
freq = fread(fileID, 'double');
fclose(fileID);
responseFlienameRegx = 'response_*.bin';
responseFilelist = dir([Directory,responseFlienameRegx]);
[~,index] = sortrows({responseFilelist.date}.');
responseFilelist = responseFilelist(index); 
clear index;

fileID = fopen([Directory, responseFilelist(1).name]);
response = fread(fileID, 'double');
fclose(fileID);
responses = zeros(length(responseFilelist), length(response));

for i = 1:length(responseFilelist)
    fileID = fopen([Directory, responseFilelist(i).name]);
    responses(i,:) = fread(fileID, 'double');
    fclose(fileID);
end

Cali_A0 = 8.615;
Cali_Aw0 = 1.180;
Cali_B0 = 9.832;
Cali_Bw0 = 1.113;

Adxs = (Cali_A0 - Axs)/Cali_Aw0*HGB0.wx(OLE.lenf);
Bdxs = (Cali_B0 - Bxs)/Cali_Bw0*HGB0.wx(OLE.lenf);
spdxs_cali = Adxs + Bdxs;
spoffs_cali = (-Adxs + Bdxs)/2;
[~,I] = min(DC_A.^2);
spoffs_cali = spoffs_cali - spoffs_cali(I); 

chirp_span  = 50;
N = length(powers);
powers = powers*chirp_span;
[DC, AC] = Match_DS(OLE, spdxs_cali, spoffs_cali);
plot(spoffs_cali, -DC, spoffs_cali, DC_A, 'o')
plot(spoffs_cali, AC, spoffs_cali, powers, 'o')

DeviceMove = -20e-6:5e-6:20e-6;
DeviceMove = [-20e-6, 20e-6];
DSMove = -20e-6:5e-6:20e-6;
data = zeros(length(DeviceMove),length(DSMove));
F0s = zeros(length(DeviceMove),length(DSMove));
F1s = zeros(length(DeviceMove),length(DSMove));
F1sdx = zeros(length(DeviceMove),length(DSMove));
for i = 1:length(DeviceMove)
    DSMove = [-20e-6:5e-6:20e-6]+DeviceMove(i);
    for j = 1:length(DSMove)
        laseroff = -DeviceMove(i);
        DSoff = DSMove(j)-DeviceMove(i);
        OLE.CalF(DSoff, mirrorGapSize);
        F0s(i,j) = OLE.F0;
        F1s(i,j) = OLE.F1;
        F1sdx(i,j) = OLE.F1*laseroff;
        data(i,j) = OLE.CoefA(laseroff, angle);
    end
end

OLE.CalF(DSoff, mirrorGapSize);
OLE.CoefA(laseroff, angle);

Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\';
fileID = fopen([Directory, 'ThCali.bin']);
ThCali = fread(fileID, 'double');
fclose(fileID);
fileID = fopen([Directory, 'ThCalifreq.bin']);
ThCalifreq = fread(fileID, 'double');
fclose(fileID);

fileID = fopen([Directory, 'fm0.bin']);
fm0 = fread(fileID, 'double');
fclose(fileID);
fileID = fopen([Directory, 'fm1.bin']);
fm1 = fread(fileID, 'double');
fclose(fileID);
fileID = fopen([Directory, 'fm2.bin']);
fm2 = fread(fileID, 'double');
fclose(fileID);
plot(1:length(fm1), fm0-fm0(1), 1:length(fm1), fm1-fm1(1), 1:length(fm1), fm2-fm2(1))

z = 0:1e-4:5e-2;
lambda = 1064*1e-9;
w0 = 60e-6;
zr = pi*w0^2/lambda;
k = 2*pi/lambda;
c1 = 1/w0^2*1./(1+(z/zr).^2);
c2 = k/2*z./(z.^2+zr^2);

Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\157-optical drive\';
FlienameRegx = 'psds_z=*.bin';
Filelist = dir([Directory,FlienameRegx]);
zs = zeros(length(Filelist), 1);
for i = 1:length(Filelist)
    r = regexp(Filelist(i).name, 'psds_z=([0-9]+\.[0-9]+).bin', 'tokens');
    zs(i) = str2num(r{1}{1});
end
[~,index] = sortrows(zs);
Filelist = Filelist(index); 
clear index;

fileID = fopen([Directory, Filelist(1).name]);
psd = fread(fileID, 'double');
fclose(fileID);li
psds = zeros(length(Filelist), length(psd));

for i = 1:length(Filelist)
    fileID = fopen([Directory, Filelist(i).name]);
    psds(i,:) = fread(fileID, 'double');
    fclose(fileID);
end

psds = psds(:, 400:600);
for i = 1:length(psds)
    plot([1:length(psds(i,:))]/10, psds(i,:),'*');
    waitforbuttonpress;
end

fm = 120000;
modeIndex = 4;
OLE.device.SetMode(fm, gamma, modeIndex)
OLE.laser.Pdc = 1e-3;
f = device.fm-50*device.gamma:0.25:device.fm+50*device.gamma;
OLE.SetInspectFreqs(f);
angles = [89.9:0.01:90.1];
angles = angles / 90*pi/2;
for angle = angles
    re = OLE.AOM_omega(angle, 20e-6);
    plot(OLE.f, 10*log10(re))
    axis([OLE.f(1) OLE.f(end) -150 0])
    waitforbuttonpress;
end

angle = 89.9/90*pi/2;
Pdcs = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
fm = 250000;
modeIndex = 8;
OLE.device.SetMode(fm, gamma, modeIndex)
f = OLE.device.fm-50*OLE.device.gamma:0.25:OLE.device.fm+50*OLE.device.gamma;
OLE.SetInspectFreqs(f);

for Pdc = Pdcs
    OLE.laser.SetPdc(Pdc)
    [re, re0] = OLE.AOM_omega(angle, 20e-6);
    plot(OLE.f, 10*log10(re),OLE.f, 10*log10(re0))
    waitforbuttonpress;
end

w0 = 120e-6;
f0 = 125e-3;
w1 = HGBeam.wz(w0, 1.064e-6, f0);
angle = (HGBeam.guoyAngle(w1, 1.064e-6, f0)+HGBeam.guoyAngle(w0, 1.064e-6, f0))/pi*180;

fm = 240000;
modeIndex = 8;
OLE.device.SetMode(fm, gamma, modeIndex)
OLE.laser.Pdc = 1e-1;
f = device.fm-50*device.gamma:0.25:device.fm+50*device.gamma;
OLE.SetInspectFreqs(f);
angles = 0:5:120;
angles = angles / 90*pi/2;
for angle = angles
re = OLE.AOM_omega(angle, 20e-6);
plot(OLE.f, 10*log10(re))
axis([OLE.f(1) OLE.f(end) -150 0])
waitforbuttonpress;
end

Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\162-chirp\';
PSDFlienameRegx = 'psds_z=*.bin';
PSDFilelist = dir([Directory,PSDFlienameRegx]);
[~,index] = sortrows({PSDFilelist.date}.');
PSDFilelist = PSDFilelist(index); 
clear index;

for i = 1:length(PSDFilelist)
    fileID = fopen([Directory, PSDFilelist(i).name]);
    psds(i,:) = fread(fileID, 'double');
    fclose(fileID);
end