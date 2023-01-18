function newOutName = updateOutname( outName, relativeOutPath, options )
% UPDATEOUTNAME renames a FloeDyn output file
% 
% FloeDyn filenames are created to be relatively unique (each one is
% appended with a randomly generated 5-digit identifier). Unfortunately
% this results in names that are not informative/reflective of model
% configurations and thus hard to parse.
%
% This function updates a given output filename from FloeDyn to be more
% informative (matching the format I've been using for input files).
% In the case of multiple files with the same name (i.e., same model
% setup), then I'll fall back on appending a random alphanumeric sequence
% to the end of the file in order to differentiate.
%
%   newOutName = updateOutname( outName )
%   newOutName = updateOutname( outName, relativeOutPath )
%   'relativeOutPath' is the file path RELATIVE to '/Floe_Cpp/io/outputs/'
%
%   newOutName = updateOutname( outName, relativeOutPath, 'Name','Value' ) 
% 
%
%   S.D.Brenner, 2022

    arguments
        outName             {mustBeText}
        relativeOutPath     {mustBeText} = '';
        options.Overwrite   (1,1) {mustBeNumericOrLogical} = 0
        options.Append      {mustBeText} = '';
        options.Prepend     {mustBeText} = 'out_';
    end

    % Load/parse data model details:
    rootDir = '~/Documents/Brown/Projects/SASIP/FloeDyn/Floe_Cpp/io/outputs/';
    fullOutPath = strcat(rootDir,relativeOutPath);
    floeOut = readFloeOut(outName,relativeOutPath);
    
    % Get model properties
    W = range( floeOut.chunk.win(1:2) );
    H = range( floeOut.chunk.win(3:4) );
    sic = floeOut.chunk.sic;
    numFloes = length(floeOut.floes.floeShapes);
    
    % Generate filename
    if log10(W)>4 && log10(H)>4
        newOutName = sprintf("%02.0fp_%ix%ikm_%if",100*sic,round(W/1e3),round(H/1e3),numFloes);
    else
        newOutName = sprintf("%02.0fp_%ix%im_%if",100*sic,round(W),round(H),numFloes);
    end

    % Optional name prepend/append
    newOutName = strcat( options.Prepend, newOutName, options.Append );

    % replace problematic characters
    newOutName = strrep(newOutName,' ','');  % remove spaces in name
    newOutName = strrep(newOutName,'.',','); % remove periods in name (replace with commas)

    % Check for repeat filename
    newOutNameXXX = newOutName;
    while isfile( strcat(fullOutPath,newOutNameXXX,'.h5') ) && ~options.Overwrite
        % generate and append random 3-digit filename addon
        charList = char([48:57,65:90,97:122]);
        XXX = charList( randi(length(charList),[1,3]) );
        newOutNameXXX = strcat(newOutName,'_',XXX);
    end
    newOutName = newOutNameXXX;

    % Rename file
    newOutName = strcat(newOutName,'.h5');   % add file extension
    if strcmp(string(outName),string(newOutName)), return; end
    [status,msg] = movefile( strcat(fullOutPath,outName), strcat(fullOutPath,newOutName) );
    if status==0; error('Error renaming file:\n%s',msg); end

end

