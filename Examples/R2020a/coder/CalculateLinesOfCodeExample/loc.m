% LOC Calculate total lines of source and generated code in a codegen run
%
% [i, o] = loc(r) 
% Takes a report information object 'r' as input, and produces two
% outputs - 'i' and 'o' containing the total lines of code in the source
% MATLAB files and generated files respectively.

% Copyright 2018 The MathWorks, Inc.
function [i, o] = loc(r)
    narginchk(1,1);

    % Assert that input is a report information object.
    assert(isa(r, 'coder.ReportInfo'), 'Input must be of type coder.ReportInfo');

    % Fetch source and generated files from the report information object.
    sourceFiles = r.InputFiles;
    generatedFiles = r.GeneratedFiles;
    
    % Count lines of code in source and generated files. Blank lines, and
    % comments are not counted.
    i = countLines(sourceFiles, true);
    o = countLines(generatedFiles, false);
end

function count = countLines(files, isSource)
    count = 0;
    for i=1:numel(files)
        f = files(i);
        if isprop(f, 'Text')
            lines = splitlines(f.Text);
            for j=1:numel(lines)
                line = strtrim(lines{j});
                if ~isempty(line) && ~isComment(line, isSource)
                    count = count + 1;
                end
            end
            clear isComment; % clear persistent variables
        end
    end
end

function result = isComment(line, isSource)
    persistent inBlockComment;
    if isempty(inBlockComment)
        inBlockComment = false;
    end
    if isSource
        result = (startsWith(line, "%") || inBlockComment);
        if line == "%{" || line == "%}"
            inBlockComment = (line ~= "%}");
        end
    else
        result = (startsWith(line, "/") || inBlockComment);
        if startsWith(line, "/*") || endsWith(line, "*/")
            inBlockComment = ~endsWith(line, "*/");
        end
    end
end