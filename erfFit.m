function re = erfFit(x, power, Aangle, tiltangle, Vmax, Mirror)
    x = x*1e-3;
    f = @(wx,P,x0,Noi,x) func(x, Aangle, tiltangle, wx, P, x0, Noi);
    [param, param_std] = fitOptions(Mirror);
    param(2) = Vmax;
    param_std(2) = 0.07*Vmax;
    param(4) = 0;
    param_std(4) = 0.1;
    re = fit(x, power', fittype(f), 'StartPoint', param, 'Lower', param - 3*param_std, 'Upper', param + 3*param_std);
    fitx = 1e-3 :1e-4 :15e-3;
    plot(fitx, re(fitx), x, power, '*');
end

function re = func(x, Aangle, tiltangle, wx, P, x0, Noi)
    re = Noi + 1/2*P*(1+erf((x-x0)*sin(Aangle)*cos(tiltangle)/(wx/sqrt(2))));
end

function [param, param_std] = fitOptions(Mirror)
    Directory = 'Z:\data\optical lever project\NORCADA_NX53515C\55-Cali\';
    CaliFlienameRegx = 'Cali*.json';
    CaliFilelist = dir([Directory,CaliFlienameRegx]);
    
    x0 = [];
    w0 = [];
    for i = 1 : length(CaliFilelist)
        fitInfo = jsondecode(fileread([Directory,CaliFilelist(i).name]));
        if fitInfo.CaliMirror == Mirror
            x0 = [x0, fitInfo.fitResult(1)];
            w0 = [w0, fitInfo.fitResult(5)];
        end
    end
    param(3) = mean(x0);
    param_std(3) = std(x0);
    param(1) = mean(w0);
    param_std(1) = std(w0);
end