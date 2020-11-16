function rdf0 = importfile_RingDown(filename, dataLines)
%IMPORTFILE2 Import data from a text file
%  RDF0 = IMPORTFILE2(FILENAME) reads data from text file FILENAME for
%  the default selection.  Returns the data as a table.
%
%  RDF0 = IMPORTFILE2(FILE, DATALINES) reads data for the specified row
%  interval(s) of text file FILENAME. Specify DATALINES as a positive
%  scalar integer or a N-by-2 array of positive scalar integers for
%  dis-contiguous row intervals.
%
%  Example:
%  rdf0 = importfile2("Z:\data\optical lever project\NORCADA_NX53515C\01-RingDown\rd_f=0.43319MHz.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 10-Nov-2020 17:16:05

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
opts.VariableNames = ["Times", "PowerdB"];
opts.VariableTypes = ["double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
rdf0 = readtable(filename, opts);

end