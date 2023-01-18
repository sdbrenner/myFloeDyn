function floeOut = floeflow( floeOut, forcing, NameValueOptions )
% FLOEFLOW calculates floE-averaged properties of the ocean floW field
%
%   floeOut = floeflow( floeOut, forcing )
%
%   The 'floeOut' structure is modified to include:
%   floeOut.ocean.com             : quantities at floe center-of-mass (CoM)
%   floeOut.ocean.avg             : floe-averaged quantities
%   floeOut.ocean.var             : variance of quanities over the floe
%   where each field is a substructure containing:
%       floeOut.ocean.(xxx).time  : elapsed time [seconds]
%       floeOut.ocean.(xxx).uo    : ocean zonal velocity
%       floeOut.ocean.(xxx).vo    : ocean meridional velocity
%       floeOut.ocean.(xxx).zeta  : ocean vertical vorticity
%       floeOut.ocean.(xxx).div   : ocean divergance
%       floeOut.ocean.(xxx).ug    : geostrophic zonal velocity
%       floeOut.ocean.(xxx).vg    : geostrophic meridional velocity
%
%   S.D.Brenner, 2022


%% Parse inputs

arguments
    floeOut                             (1,1) {mustBeA(floeOut,"struct")}
    forcing                             (1,1) {mustBeA(forcing,"struct")}
    NameValueOptions.outputFloeMasks    (1,1) {mustBeNumericOrLogical} = false;
    NameValueOptions.outputDomainMask   (1,1) {mustBeNumericOrLogical} = false;
end


%% Subset floeOut data to put on common timestep with forcing
% (assumes lower temporal resolution for ocean forcing)

tol = 120; % timestep tolerance [seconds]
[LIA,~] = ismembertol( floeOut.state.time, forcing.t,tol,'DataScale',1 );
floeOut = subsetFloeOut( floeOut, LIA );

%% Create binary masks 
% Binary masks are used for averaging the forcing fields over each floe. 
%
% I SHOULD return these as part of the 'floeOut' struct; however, even as
% sparse matrices they significantly bloat the filesize of 'floeOut'. In
% order to keep 'floeOut' smaller (for saving/reloading for batch-script
% analysis), the default behaviour will be to not output the masks, but
% they're included on request.

if ~isfield(  floeOut.floes,'floeMasks' )
%     disp('Making binary masks for ice floes...');
    floeOutMask = makeFloeMask( floeOut, forcing.x,forcing.y );
else
    floeOutMask = floeOut;
end
[M,L] = size( floeOutMask.floes.floeMasks );

% Save floe masks in struct if requested
if NameValueOptions.outputFloeMasks
    floeOut = floeOutMask;
end

% Add domain floe mask if requested
% (this is basically "free" once the floeMasks are calculated)
if NameValueOptions.outputDomainMask
    nx = length(forcing.x);
    ny = length(forcing.y);
    for l = 1:L
        floeMasksL = floeOutMask.floes.floeMasks(:,l);
        fml = reshape( full( [floeMasksL{:}] ), nx,ny,M);
        dml = sparse(any(fml,3));
        floeOut.chunk.iceMask{l} = dml;
    end
end

%% Calculate/interpolate ocean forcing properties
disp('Calc. ocean properties');

forcing = renameStructField(forcing,'u','uo');
forcing = renameStructField(forcing,'v','vo');

nx = length(forcing.x); ny = length(forcing.y);
geoTrue = isfield(forcing,'ug');
% Define fields for averaging
frcFlds = ["uo","vo","zeta","div"];
if geoTrue
    frcFlds(end+1:end+2) = ["ug","vg"];
end

t = ones(M,1)*floeOut.state.time;
xp = floeOut.state.xp;
yp = floeOut.state.yp;

numFlds = length(frcFlds);
% Loop through forcing fields
for n = 1:numFlds
    % Floe Center-of-mass values
    floeOut.ocean.com.(frcFlds(n)) = reshape( interp3( forcing.x, forcing.y, forcing.t,...
                                        forcing.(frcFlds(n)), xp(:),yp(:),t(:) )  ,M,L);
%     floeOut.ocean.com.(frcFlds(n)) = reshape( floeOut.ocean.com.(frcFlds(n)),M,L);
    % Loop through time
    floeOut.ocean.com.(frcFlds(n)) = NaN(M,L);
    floeOut.ocean.avg.(frcFlds(n)) = NaN(M,L);
    floeOut.ocean.var.(frcFlds(n)) = NaN(M,L);
    for l = 1:L
        k = l + (n-1)*L;
%         fprintf('floeflow calc-step %3g of %3g...\n',k,numFlds*L);     
%         % Floe Center-of-mass values
%         floeOut.ocean.com.(frcFlds(n))(:,l) = interp2( forcing.x, forcing.y,...
%             forcing.(frcFlds(n))(:,:,l), floeOut.state.xp(:,l), floeOut.state.yp(:,l) );
        % Floe-averaged quantities & variance
        floeMasksL = floeOutMask.floes.floeMasks(:,l);
        floeMasksL = reshape( ([floeMasksL{:}]), nx*ny,M );
%         % Average all floes at once:
%         nzv = -100; %ensure valML is always nonzero
%         valML = floeMasksL.*reshape( forcing.(frcFlds(n))(:,:,l)+nzv,nx*ny, [] ); % [this operation is SLOW]
%         mInds = nonzeros(floeMasksL .* (1:M));
%         floeOut.ocean.avg.(frcFlds(n))(:,l) = accumarray( mInds, nonzeros(valML)-nzv,[],@mean  );
%         floeOut.ocean.var.(frcFlds(n))(:,l) = accumarray( mInds, nonzeros(valML)-nzv,[],@var  );
        % Loop through floes and average
        valML = reshape( forcing.(frcFlds(n))(:,:,l), nx*ny, [] );
        for m = 1:M
            floeOut.ocean.avg.(frcFlds(n))(m,l) = mean( valML(floeMasksL(:,m)) );
            floeOut.ocean.var.(frcFlds(n))(m,l) =  var( valML(floeMasksL(:,m)) );
        end
    end
end
    

end

%% Calculate/interpolate ocean forcing properties

%{
geoTrue = isfield(forcing,'ug');


floeOut.ocean.time = floeOut.state.time; 
% Loop through time
for l = 1:L
    floeX = floeOut.state.xp(:,l);
    floeY = floeOut.state.yp(:,l);
    fuL = forcing.u(:,:,l);
    fvL = forcing.v(:,:,l);
    fzL = forcing.zeta(:,:,l);
    fdL = forcing.div(:,:,l);
    % 2D interpolate to find point values
    floeOut.ocean.com.uo(:,l) = interp2( forcing.x, forcing.y, fuL, floeX, floeY );
    floeOut.ocean.com.vo(:,l) = interp2( forcing.x, forcing.y, fvL, floeX, floeY );
    floeOut.ocean.com.zeta(:,l) = interp2( forcing.x, forcing.y, fzL, floeX, floeY );
    floeOut.ocean.com.div(:,l) = interp2( forcing.x, forcing.y, fdL, floeX, floeY );
    % include geostrophic velocities if they're in the file
    if geoTrue
        fugL = forcing.u(:,:,l);
        fvgL = forcing.v(:,:,l);
        floeOut.ocean.com.uo(:,l) = interp2( forcing.x, forcing.y, fugL, floeX, floeY );
        floeOut.ocean.com.vo(:,l) = interp2( forcing.x, forcing.y, fugL, floeX, floeY );
    end

    % Loop through floes for average/variance
    for m = 1:M
        % averages
        floeMaskML = floeOutMask.floes.floeMasks{m,l};
        floeOut.ocean.avg.uo(m,l) = mean( fuL(floeMaskML) );
        floeOut.ocean.avg.vo(m,l) = mean( fvL(floeMaskML) );
        floeOut.ocean.avg.zeta(m,l) = mean( fzL(floeMaskML) );
        floeOut.ocean.avg.div(m,l) = mean( fdL(floeMaskML) );
        % variance
        floeOut.ocean.var.uo(m,l) = var( fuL(floeMaskML) );
        floeOut.ocean.var.vo(m,l) = var( fvL(floeMaskML) );
        floeOut.ocean.var.zeta(m,l) = var( fzL(floeMaskML) );
        floeOut.ocean.var.div(m,l) = var( fdL(floeMaskML) );
        if geoTrue % include geostrophic velocities if they're in the file
            % average
            floeOut.ocean.avg.ug(m,l) = mean( fugL(floeMaskML) );
            floeOut.ocean.avg.vg(m,l) = mean( fvgL(floeMaskML) );
            % variance
            floeOut.ocean.var.ug(m,l) = var( fugL(floeMaskML) );
            floeOut.ocean.var.vg(m,l) = var( fvgL(floeMaskML) );
        end
    end
end

end
%}




% 
% function [meanX, varX] = sparseStat( X, mask )
%     meanX = sum(X)./sum(mask);
%     mXS = mask.*meanX;
%     Y = X;
%     Y(mask) = (X(mask)-mXS(mask)).^2;
%     varX = sum(Y)./(sum(mask)-1);
% end