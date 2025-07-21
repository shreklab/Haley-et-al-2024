function [] = writeSourceData(fileName,header,tableNames,dataTables,notes)
% [] = WRITESOURCEDATA(FILENAME,HEADER,TABLENAMES,DATATABLES,NOTES)
%
%   WRITESOURCEDATA compiles and writes source data, including a header,
%   multiple data tables with titles, and accompanying notes, into a
%   single, formatted Excel (.xlsx) file. The function sequentially appends
%   each table and its metadata, making it suitable for creating
%   comprehensive source data files for figures or analyses.
%
%   INPUTS:
%       - fileName [char]: The full path and file name for the output
%           .xlsx file.
%       - header [cell]: A cell array of strings containing the header
%           information to be written at the top of the file.
%       - tableNames [cell]: A cell array of strings, where each string
%           serves as the title for a corresponding data table.
%       - dataTables [cell]: A cell array where each element is a MATLAB
%           'table' containing the data to be written. The order must
%           correspond to tableNames.
%       - notes [cell]: A cell array where each element is a cell array of
%           strings containing notes to be written to the right of the
%           corresponding data table.
%
%   OUTPUTS:
%       This function does not return any variables to the MATLAB
%       workspace. Its primary output is the .xlsx file saved to the path
%       specified by fileName.
%
%   Written 7/11/2025 by Jess Haley in MATLAB R2024a.
%
%   See also WRITECELL, WRITETABLE, READMATRIX.

% Write deader
writecell(header,fileName,'AutoFitWidth',false,...
    'Range',[char(64+1),num2str(1)],'WriteMode','replacefile')

for i = 1:length(tableNames)
    % Check file size
    rowNum = size(readmatrix(fileName,'OutputType','char','Range','A1'),1);

    % Write table name
    writecell(tableNames(i),fileName,'AutoFitWidth',false,...
        'Range',[char(64+1),num2str(2+rowNum(1))]);

    % Write table
    writetable(dataTables{i},fileName,'AutoFitWidth',false,...
        'Range',[char(64+2),num2str(4+rowNum(1))],'WriteRowNames',true)

    % Write notes
    if strcmp(dataTables{i}.Properties.DimensionNames{1},'Row')
        offset = 3;
    else
        offset = 4;
    end
    if ~isempty(notes{i})
        writecell(notes{i},fileName,'AutoFitWidth',false,...
            'Range',[char(64+offset+width(dataTables{i})),num2str(4+rowNum(1))])
    end
end

end