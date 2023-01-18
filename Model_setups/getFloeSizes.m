function floeSizes = getFloeSizes( FSD,numFloes,totalIceArea )
% GETFLOESIZES generates a list of floe sizes for a given FSD satisfying
% either a total number of floes or a total sea ice area
%
%   floeSizes= getFloeSizes( FSD,numFloes,totalIceArea )
%   [floeSizes,ffun,Ffun] = getFloeSizes( FSD,numFloes,totalIceArea )
%
%   S.D.Brenner, 2022

    %% Parse inputs
    % add error checking

    % Check if upper and lower bounds of FSD are equal
    if FSD(2)==FSD(1) 
        floeSizes = getUniformFloeSizes( FSD,numFloes,totalIceArea );
%         ffun = @(x) NaN*x; Ffun = @(x) NaN*x;
        return;
    end


    %% Get FSD
    % Floe size distribution is specified based on a truncated power-law,
    % following Stern et al., 2018 (https://doi.org/10.1525/elementa.304)
    % so that cumulative distribution is concave down

    % Get FSD parameters from input
    a = FSD(1);     % mininum floe size
    b = FSD(2);     % maximum floe size
    alf = FSD(3);   % non-cumulative power-law exponent
    % parameters associated with distribution (see Stern2018)
    B = (b/a).^(1-alf);
    c = (alf-1)*a^(alf-1)/(1-B);
    C = c/(alf-1);
    R = C*b^(1-alf);
    % Theoretical FSD functions (see Stern2018)
%     ffun = @(x) c*x.^(-alf);
%     Ffun = @(x) C*x.^(1-alf)-R;
%     [ffun,Ffun] = fsdFuns( a,b,alf ); 

    %% Determine number of floes (if not specified)

    if isempty(numFloes) && ~isempty(totalIceArea)
        % Determine number of floes necessary:
        % [ Create a curve of number of floes vs sea ice area , then
        %   numerically invert to find the number of floes necessary
        %   for the requested total sea ice area ]
        maxN = 5e3;
        nF = 0:100:maxN;
        M = numel(nF);
        iceArea = zeros(1,M);
        for m = 1:M
            k = linspace(0,nF(m),nF(m));    
            Fdis = 1-k/nF(m);
            xdis = ((Fdis+R)/C).^(1/(1-alf));
            iceArea(m) = sum( xdis.^2 );
        end
        numFloes = round( interp1( iceArea, nF, totalIceArea,"pchip",'extrap' ) );
    end

    %% Generate floes

    k = linspace(0,numFloes,numFloes);    
    Fdis = 1-k/numFloes;
    xdis = ((Fdis+R)/C).^(1/(1-alf));

    %% Re-adjust sizes (if target ice area specified)

    if nargin>2 && ~isempty(totalIceArea)
        xdis = xdis*sqrt( totalIceArea/sum(xdis.^2) ); 
    end

    %% Organize output 

    floeSizes = sort(xdis,"descend");

end













function floeSizes = getUniformFloeSizes( FSD,numFloes,totalIceArea )

    floeSize = FSD(1);

    % Determine number of floes (if not specified) 
    if isempty(numFloes) && ~isempty(totalIceArea)
        numFloes = round( totalIceArea/floeSize^2 );
    end

    % Create vector containing repeated values of floeSize:
    floeSizes = repmat(floeSize,[1,numFloes]);

end

