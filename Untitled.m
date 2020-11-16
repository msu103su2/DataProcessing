folder = 'D:/temp/';
fileNameRegex = 'rd_f=([0-9].[0-9]+)MHz.csv';
filelist = dir([folder, '*.csv']);
FQ = zeros(size(filelist,1),2);
for i = 1:size(filelist,1)
    if regexp(filelist(i).name, fileNameRegex)
        tokens = regexp(filelist(1).name, fileNameRegex, 'tokens');
        FQ(i,1) = str2double(tokens{1}{1});
        preRD = FQ(1:CloestIndex(FQ(:,1),1.5),2);
    end
end
function Idx = CloestIndex(array,element)
    [~,Idx] = min(abs(array - element));
end

dis