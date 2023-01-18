function [kr,Zr] = shellSpectra2D( x,y,z )
% 
%   [k,Y] = shellSpectra2D( x,y,Z )
%
%   S.D.Brenner, 2022

    nx = length(x);
    ny = length(y);
    dx = range(x)/(nx-1);
    dy = range(y)/(ny-1);
    
    % Get power
    Z = fft2( z ); 
    Z = fftshift(Z);
    Z = 2*(dx/nx)*(dy/ny)*abs(Z).^2;
    
    % Get wavenumbers
    kx = linspace(-1/(2*dx),1/(2*dx),nx);
    ky = linspace(-1/(2*dy),1/(2*dy),ny);
    dkx = (1/dx)/(nx-1);
    dky = (1/dy)/(ny-1);
    dk = sqrt( dkx.^2 + dky.^2 );
    [kX,kY] = meshgrid(kx,ky);
    [~,kR] = cart2pol(kX,kY);
    
    % Average in "shells"
    % [see e.g., XXXX ]  
%     krN = min( [ 1/(2*dx), 1/(2*dy)] ); % (sort-of Nyquist wavenumber?)
%     kBin = 1*dk;
%     bins = (dk/2):kBin:( krN+dk );
    krN = min( [ 1/(3*dx), 1/(3*dy)] ); % (sort-of Nyquist wavenumber?)
    kBin = 1*dk;
    bins = (0):kBin:( krN+dk );   
    [kr,Zr] = binFun( kR(:),Z(:),bins,'returnX','binCenters'); % radially average
    kr = kr(:); Zr = Zr(:);

    % Rescale? (spectral density)
    Zr = (pi/2) * kr .* Zr;



end


