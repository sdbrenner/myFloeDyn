%% Create HD5 forcing file
% Check to see if filesize is notably smaller
% ...then try to figure out if FloeDyn read-in is any faster

%% Clean workspace

clearvars -except forcing;
clc;
close all;

%% Add paths

addpath('../Process_scripts/');

%% Load forcing data

forcePath = 'forcings/';
forceName = '_input_forcing_SMgeo_5day.mat';
forcing = readForcing(forceName,forcePath);

forcing = rmfield(forcing,'div','zeta');

% save('fTest.mat','forcing','-struct','-v7');

%% Create HD5 file

fName = 'fh5Test.h5';
if exist(fName,'file'), delete(fName); end

% Create empty file
flds = fields(forcing);
for n = 1:7
    sz = size( forcing.(flds{n}) );
    h5create( fName, ['/',flds{n}], sz,'Datatype','single','Chunksize',sz,'Deflate',9);
end

% Write data
for n = 1:7
    sz = size( forcing.(flds{n}) );
    h5write( fName, ['/',flds{n}], single(forcing.(flds{n})) );
end

% h5fo = h5info(fName);

%% 
% Try the struct2hdf5 from the FEX
% (https://www.mathworks.com/matlabcentral/fileexchange/106470-struct2hdf5)


% struct2hdf5(forcing,'forcing','','fh5Test.h5');



function struct2hdf5(data,dataset,fp,fn)
% Author: Isaac Li
% This code is based on sample code from: 
% https://support.hdfgroup.org/HDF5/examples/api18-m.html
% Version 20220209:
% The function takes a struct and write to a compound HDF5 file. the fields
% in the struct will automatically populate in the HDF5.
% Automatically handles single, double, uint 8,16,32 data types, feel free
% to add more in the case structure below if your data has more.
% 20220715 - bug fix
% data - struct
% dataset - dataset name
% fp - filepath
% fn - filename
% This file is intended for use with HDF5 Library version 1.8
	full_fn = [fp '\' fn];
% 	dims    = length(data.frame);
	fn		= fieldnames(data);
	dims	= length(getfield(data,fn{1}));
	%% Create a new file using the default properties.
	file = H5F.create (full_fn, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
	%% Create the required data types for each field
	data_fieldnames = fieldnames(data);
	Nfield = length(data_fieldnames);
	for k = 1:Nfield
		switch class(data.(data_fieldnames{k}))
			% For more data types, see:
			% https://support.hdfgroup.org/HDF5/doc1.8/RM/PredefDTypes.html
			case 'single'
				h5_datatype(k) = H5T.copy('H5T_NATIVE_FLOAT');
			case 'double'
				h5_datatype(k) = H5T.copy('H5T_NATIVE_DOUBLE');
			case 'uint32'
				h5_datatype(k) = H5T.copy('H5T_NATIVE_UINT32');
			case 'uint16'
				h5_datatype(k) = H5T.copy('H5T_NATIVE_UINT16');
			case 'uint8'
				h5_datatype(k) = H5T.copy('H5T_NATIVE_UINT8');
			otherwise
				error('struct2hdf5: datatype not recognized, add when needed');
		end
		sz(k) = H5T.get_size(h5_datatype(k));
	end
	offset(1)=0;
	offset(2:Nfield) = cumsum(sz(1:(Nfield-1)));
	%% Create the compound datatype for memory.
		memtype = H5T.create ('H5T_COMPOUND', sum(sz));
		for k = 1:Nfield
			H5T.insert (memtype,data_fieldnames{k},offset(k),h5_datatype(k));	
		end
	%% Create the compound datatype for the file.  Because the standard
	% types we are using for the file may have different sizes than
	% the corresponding native types, we must manually calculate the
	% offset of each member.
		filetype = H5T.create ('H5T_COMPOUND', sum(sz));
		for k = 1:Nfield
			H5T.insert (filetype,data_fieldnames{k},offset(k),h5_datatype(k));	
		end
	%% Create dataspace.  Setting maximum size to [] sets the maximum
	% size to be the current size.
		space = H5S.create_simple (1,fliplr(dims), []);
	%% Create the dataset and write the compound data to it.
		dset = H5D.create (file, dataset, filetype, space, 'H5P_DEFAULT');
		H5D.write (dset, memtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', data);
	%% Close and release resources.
		H5D.close(dset);
		H5S.close(space);
		H5T.close(filetype);
		H5F.close(file);
	%% debug
	% 	h5disp(full_fn);
end
