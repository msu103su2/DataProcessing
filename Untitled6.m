Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\82-chirp_over_dy\';

PSDFlienameRegx = 'PSD_*.bin';
PSDFilelist = dir([Directory,PSDFlienameRegx]);
[~,index] = sortrows({PSDFilelist.date}.');
PSDFilelist = PSDFilelist(index); 
clear index;

DrivePSDFlienameRegx = 'DrivePSD_*.bin';
DrivePSDFilelist = dir([Directory,DrivePSDFlienameRegx]);
[~,index] = sortrows({DrivePSDFilelist.date}.');
DrivePSDFilelist = DrivePSDFilelist(index); 
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
    fileID = fopen([DrivePSDFilelist(i).folder,'\',DrivePSDFilelist(i).name]);
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

rho = 3100;
a = 1.5e-3;
b = 3.5e-3;
h = 100e-9;
w0 = 60e-6;
flen = 600e-3;
lambda = 1064e-9;
sp = SplitPair(0, 0);
device = Membrane(a, b, h, rho);
device.SetMode(182194, 5, 2);
laser = Laser(lambda);
laser.SetPdc(1.7e-3);
laser.SetAC(1.02e-3, 0.1);
PD = PhotoDetector(300);
OLE = OLExperiment(sp, device, laser, w0, flen, PD);
angle = HGBeam.guoyAngle(w0, lambda, flen);

f = device.fm-500*device.gamma:0.1*device.gamma:device.fm+500*device.gamma;
OLE.SetInspectFreqs(f);
OLE.CalF(0, 0);

spdxs = 0:0.01e-5:2e-5;
distance = zeros(1,length(spdxs));
[total, signal, noise] = OLE.GetPSDs(1e-5, angle);
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

Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\114-LineScan\';
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

Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\135-Dtest\';
fileID = fopen([Directory, 'DC_A.bin']);
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