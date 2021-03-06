function [DC, AC] = Match_DS(OLE, spdxs, spoffs)
    F0s = zeros(1, length(spdxs));
    F1s = zeros(1, length(spdxs));
    for i = 1:length(spdxs)
        OLE.CalF(spoffs(i), spdxs(i));
        F0s(i) = OLE.F0;
        F1s(i) = OLE.F1;
    end
    DC = (OLE.laser.Pdc*F0s*OLE.PD.WtoA*OLE.PD.Gain*1e3);
    AC = (OLE.laser.Pac*F0s*OLE.PD.WtoA*OLE.PD.Gain*1e3).^2/2;
end