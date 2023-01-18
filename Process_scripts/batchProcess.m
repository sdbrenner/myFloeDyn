% Batch process 

%% Clean workspace
clear -except forcing;
clc;
close all;

%% Add paths

addpath('../Analysis/');


%% Load forcing details

forceName = '_input_forcing_SMgeo_5day.mat';
% forceName = '_input_forcing_2DQGn.mat';
forcePath = 'forcings/';
if ~exist('forcing','var')
    forcing = readForcing(forceName,forcePath);
end

%% Update filenames

fDir = 'sic_batches/FSD_300_7500_2_geo/';
% fDir = 'sic_batches/FSD_3000_30000_2_geo/';
files = dir(['../../Floe_Cpp/io/outputs/',fDir,'*.h5']);
if isempty(files)
    error('No files; check directory name.');
end
% sort by newest-first
fTime = datenum({files.date});
[~,ind] = sort(fTime,'descend');
files = files(ind);
% loop and rename (saving list of orginal names for reference)
numFiles = length(files);
nameTable = strings(numFiles,2);
for n = 1:numFiles
%     appStr = datestr(now,'mmdd');
%     appStr = datestr(files(n).date,'mmdd');
    nameTable(n,1) = string(files(n).name); 
    nameTable(n,2) = updateOutname(files(n).name, fDir, Overwrite=true );
end

%% Loop through files in directory

files = dir(['../../Floe_Cpp/io/outputs/',fDir,'*.h5']);
numFiles = length(files);
for n = numFiles:-1:1 % go backwards to allocate structure memory
    fName = strcat(fDir,files(n).name);
    fprintf('Processing %s...\n',fName)
    oldNameN = nameTable( strcmp( files(n).name, nameTable(:,2) ), 1 );
    % read and process
    disp('Read file...');
    fO = readFloeOut( files(n).name,fDir,...
                      ...orginalOutName = oldNameN,...
                      forceName = strcat(forcePath,forceName) );
    fO = floeflow( fO, forcing,outputFloeMasks=true,outputDomainMask=true ); % include binary masks for IOBL calcs
    disp('Calcuate IOBL dissipation...');
    fO = calcIOBLDiss( fO, forcing );
    % remove binary masks to reduce variable sizes
    fO.floes = rmfield(fO.floes,'floeMasks'); 
    fO.chunk = rmfield(fO.chunk,'iceMask');
    % add to master structure
    floeOut(n) = fO;
    clc;
end


%% Save floeOut output

% Generate filename
sInd = regexp(fDir(1:end-1),'/');
if ~isempty(sInd)
    sInd = sInd(end)+1; 
else
    sInd = 1;
end
sName = strcat('floeOut_',fDir(sInd:end-1) );

% Save
save( strcat('../MAT_Files/',sName),'floeOut');
fprintf('Saved %s\n',strcat('../MAT_Files/',sName));