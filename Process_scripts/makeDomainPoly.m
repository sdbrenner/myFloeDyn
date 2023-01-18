function floeOut = makeDomainPoly( floeOut )
% MAKEDOMAINPOLY returns the MATLAB polyshapes for the assembly of floes
%
%   floeOut = makeDomainPoly( floeOut ), where 'floeOut' was
%   produced by the 'readFloeOut' function
%   The 'floeOut' structure is modified to include a new field:
%       floeOut.chunk.icePoly
%  
%   The polygon accounts for periodic BCs.
%
%   S.D.Brenner, 2022

% Check if the floeXYT field already exists; if not, run calcFloeXYT
if ~ isfield( floeOut.floes, 'floeXYT' )
    floeOut = calcFloeXYT( floeOut );
end


% Parse domain size from forcing floeDyn
X = floeOut.chunk.win(1);
Y = floeOut.chunk.win(3);
W = range( floeOut.chunk.win(1:2) );
H = range( floeOut.chunk.win(3:4) );

% Extract data
floeXYT = floeOut.floes.floeXYT;
[M,L] = size(floeXYT);

% Periodic BCs
if floeOut.state.PBC
    % Create x/y offset for all possible translations
    % (in normalized nx/ny distances)
    xOff = W*reshape(repmat(-1:1,3,1) ,1,[]);
    yOff = H*reshape(repmat(-1:1,3,1)',1,[]);
else
    xOff = 0;
    yOff = 0;
end

% Supress polyshape warning
wIdentif = 'MATLAB:polyshape:repairedBySimplify';
warning('off',wIdentif);

% Loop through time
icePoly = cell(1,L);
for l = 1:L
    allVert = [];
    % Loop through floes
    for m = 1:M
        % Get vertices for floe polygon (+ periodic multiples)
        vnX = [ floeXYT{m,l}(1,:).' ; NaN ] + xOff ;
        vnY = [ floeXYT{m,l}(2,:).' ; NaN ] + yOff ;
        % Only retain floe multiples if they contain a vertex point in domain
        cs = any( vnX>=X & vnX<=(X+W) & vnY>=Y & vnY<=(Y+H) ,1);
        vnX = vnX(:,cs); vnY = vnY(:,cs);
        % Create combined [X,Y] vertex matrix 
        vertML = [vnX(:),vnY(:)];
        allVert = [allVert; vertML];
    end
    % Create combined polyshape for whole domain
    icePoly{l} = polyshape(allVert,KeepCollinearPoints=true,Simplify=true);
end

% Add to structure
floeOut.chunk.icePoly = icePoly;

end


%% DEPRICATED CODE
% Older, slower loop; I'm keeping it copied here just in case 

% % Loop through time 
% floeMasks = cell(M,L);
% for l = 1:L
%     % Loop through floes 
%     for m = 1:M
%         % Create mask (if PBC, need to add "wraparound")
%         % [ there's probably a better way to do this ]
%         if floeOut.state.PBC 
%             px = ceil( (max( floeSize(m) )/W)*nx );
%             py = ceil( (max( floeSize(m) )/H)*ny );
%             vertNormMK = (floeXYT{m,l}.'-[X,Y])./[dx,dy]+[1,1]+[px,py];
%             imL = poly2mask(vertNormMK(:,1),vertNormMK(:,2),nx+2*px,ny+2*py);
%             % add x-wraparound
%             imL( :, px+(1:px) ) = imL( :, px+(1:px) ) + imL( :, px+nx+(1:px) ); % R2L
%             imL( :, nx+(1:px) ) = imL( :, nx+(1:px) ) + imL( :, (1:px) );       % L2R
%             % add y-wraparound
%             imL( py+(1:py),: ) = imL( py+(1:py),: ) + imL( py+ny+(1:py),: ); 
%             imL( ny+(1:py),: ) = imL( ny+(1:py),: ) + imL( (1:py),: );       
%             % extral central mask
%             mask = imL( px+(1:nx) , py+(1:ny) );
%         else
%             vertNormMK = (floeXYT{m,l}.'-[X,Y])./[dx,dy]+[1,1];
%             mask = poly2mask(vertNormMK(:,1),vertNormMK(:,2),nx,ny);
%         end
%         floeMasks{m,l} = sparse(mask);
%     end
% end   
