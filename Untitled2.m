folder = 'Z:/data/optical lever project/Die 000_13/04-seachModesForAll_withFineOnModes/';
flienameRegx = '01CF*.csv'
filelist = dir([folder,flienameRegx]);
Qs = zeros(size(filelist,1),1);
Ws = zeros(size(filelist,1),1);
Fs = zeros(size(filelist,1),1);
for i = 1 : size(filelist,1)
    newTable = importfile1([folder,filelist(i).name]);
    [Q, fmax, Lindex, Rindex] = cq(newTable.FrequenciesHz, newTable.PSDdB2Hz);
    Qs(i) = Q;
    Fs(i) = fmax;
    Ws(i) = newTable.FrequenciesHz(Rindex) - newTable.FrequenciesHz(Lindex);
    half_F = [newTable.FrequenciesHz(Lindex),newTable.FrequenciesHz(Rindex)];
    half_P = [newTable.PSDdB2Hz(Lindex),newTable.PSDdB2Hz(Rindex)];
    plot(newTable.FrequenciesHz, newTable.PSDdB2Hz, half_F, half_P,'*');
    legend(['Q = ', num2str(Q)],['width = ', num2str(half_F(2) - half_F(1)),'@',num2str(Q*(half_F(2) - half_F(1)))]);
end
figure('Name','new2');
index = Qs > 50e3;
Qs = Qs(index);
Fs = Fs(index);

plot(1, Fs,'*');

function [q,fmax, Lindex,Rindex] = cq(f,p)
    [pmax,pmaxI] = max(p);
    pHalfMax = pmax - 3;
    fmax = f(pmaxI);
    [temp, Lindex] = min(abs(p(1:pmaxI)-pHalfMax));
    fLeft = f(Lindex);
    [temp, Rindex] = min(abs(p(pmaxI:end)-pHalfMax));
    Rindex = Rindex + pmaxI - 1;
    fRight = f(Rindex);
    FWHM = fRight - fLeft;
    q = fmax/FWHM;
end