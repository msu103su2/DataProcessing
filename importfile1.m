function test351 = importfile1(filename, dataLines)
%IMPORTFILE1 Import data from a text file
%  TEST351 = IMPORTFILE1(FILENAME) reads data from text file FILENAME
%  for the default selection.  Returns the data as a table.
%
%  TEST351 = IMPORTFILE1(FILE, DATALINES) reads data for the specified
%  row interval(s) of text file FILENAME. Specify DATALINES as a
%  positive scalar integer or a N-by-2 array of positive scalar integers
%  for dis-contiguous row intervals.
%
%  Example:
%  test351 = importfile1("Z:\data\optical lever project\pull\000_13_test_35_1.019~1.022MHz.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 29-Jun-2020 11:34:25

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["FrequenciesHz", "PSDdB2Hz"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
test351 = readtable(filename, opts);

end