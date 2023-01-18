function makeInput(pShape,window,inName)


    %% Loop through poly matrix and generate file inputs
    
    numFloes = length(pShape);

    for k = 1:numFloes

        [fsXY(1),fsXY(2)] = centroid(pShape(k));
    
        % Create variables used for making the HDF5 file
        floeList{k} = sprintf('%u',k-1);
        floeXY{k} = flipud((pShape(k).Vertices-fsXY)).' ;
        floe_states(1:2,k) = fsXY;
        floe_states(8:9,k) = fsXY;

    end




    %% Create h5 file:

    % if the file extension isn't already given, add it:
    inName = convertStringsToChars(inName);
    if ~strcmp(inName(end-2:end),'.h5')
        inName = strcat(inName,".h5");   % add file extension 
    end


    fPath = "../Floe_Cpp/io/inputs/custom/";
    fName = strcat(fPath,inName);
    if exist(fName,'file'), delete(fName); end
    
    % Create empty HDF5 file
    h5create(fName,'/floe_states',size(floe_states,1:3),'ChunkSize',size(floe_states,1:3));
    h5create(fName,'/window',length(window) );
    for k = 1:length(floeList)
    %     S = hfo.Groups.Datasets(k).Dataspace.Size;
        S = size(floeXY{k});
        h5create(fName, ['/floe_shapes/',floeList{k}], S );
    end
    % Write data to file
    h5write(fName,'/floe_states',floe_states);
    h5write(fName,'/window',window);
    for k = 1:length(floeList)
        % check overwrite size
    %     S = hfo.Groups.Datasets(k).Dataspace.Size(2);
        h5write(fName, ['/floe_shapes/',floeList{k}], floeXY{k} );
    end
    fprintf('Saved File: %s\n',fName)



end