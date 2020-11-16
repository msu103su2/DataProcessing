folder = 'Z:\data\optical lever project\NORCADA_NX53515C\01-RingDown\';
fileNameRegex = 'rd_f=([0-9].[0-9]+)MHz.csv';
filelist = dir([folder, '*.csv']);
FQ = zeros(size(filelist,1),4);
for i = 1:size(filelist,1)
    if regexp(filelist(i).name, fileNameRegex)
        tokens = regexp(filelist(i).name, fileNameRegex, 'tokens');
        FQ(i,1) = str2double(tokens{1}{1})*1e6;
        
        data = importfile_RingDown([folder, filelist(i).name]);
        preRDEndIdx = CloestIndex(data.Times,1.5);
        preRD = data.PowerdB(1:preRDEndIdx);
        preRD_mean = mean(preRD);
        preRD_std = std(preRD);
        
        lineStart = CloestIndex(data.Times,1.85);
        for j = CloestIndex(data.Times,1.85) : size(data,1)
            if data.PowerdB(j) < preRD_mean - 3*preRD_std
                lineStart = j;
                break;
            end
        end
        
        lineEnd = lineStart;
        for j = lineStart : size(data,1)
            if data.PowerdB(j) < -80
                lineEnd = j;
                break;
            end
            if j == size(data,1)
                lineEnd = j;
            end
        end
        
        dataToFit = data(lineStart:lineEnd,:);
        x = dataToFit.Times;
        y = dataToFit.PowerdB;
        p = polyfit(x, y, 1);
        FQ(i,2) =-2*pi*FQ(i,1)/p(1);
        FQ(i,4) = p(1);

        yfit = polyval(p, x);
        SSresid = sum((y - yfit).^2);
        SStotal = (length(y)-1)*var(y);
        rsq = 1 - SSresid/SStotal;
        FQ(i,3) = rsq;
        %{
        tiledlayout(2,1)
        nexttile;
        plot(data.Times,data.PowerdB);
        nexttile;
        plot(x,y,'*',x,yfit,'-k')
        xlabel({'f=',num2str(FQ(i,1))});
        waitforbuttonpress;
        %}
    end
end
function Idx = CloestIndex(array,element)
    [~,Idx] = min(abs(array - element));
end
