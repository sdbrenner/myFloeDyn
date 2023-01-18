function floeOut = makeFloeMask( floeOut, domainX,domainY )
% MAKEFLOEMASK returns the binary masks corresponding to each of the floes
%
%   floeOut = makeFloeMask( floeOut, domainX,domainY ), where 'floeOut' was
%   produced by the 'readFloeOut' function, and 'domainX' and 'domainY' are
%   each vectors specifying the x,y coordinates of the domain that the
%   masks are to be created in (e.g, domainX = forcing.x ).
%   The 'floeOut' structure is modified to include a new field:
%       floeOut.floes.floeMask
%  
%   Floe masks account for periodic BCs.
%
%   S.D.Brenner, 2022

    % Check if the floeXYT field already exists; if not, run calcFloeXYT
    if ~ isfield( floeOut.floes, 'floeXYT' )
        floeOut = calcFloeXYT( floeOut );
    end
    
    % Parse domain size/grid data from forcing domain
    Xf = min(domainX);
    Yf = min(domainY);
    Wf = range(domainX);
    Hf = range(domainY);
    nx = length(domainX);
    ny = length(domainY);
    dx = Wf/(nx-1);
    dy = Hf/(ny-1);
    % Parse domain size from forcing floeDyn
    % Xi = floeOut.chunk.win(1);
    % Yi = floeOut.chunk.win(3);
    Wi = range( floeOut.chunk.win(1:2) );
    Hi = range( floeOut.chunk.win(3:4) );
    
    
    % Extract data
    floeXYT = floeOut.floes.floeXYT;
    [M,L] = size(floeXYT);
    
    % Periodic BCs
    if floeOut.state.PBC
        % Create x/y offset for all possible translations
        % (in normalized nx/ny distances)
        xOff = (Wi/dx)*reshape(repmat(-1:1,3,1) ,1,[]);
        yOff = (Hi/dy)*reshape(repmat(-1:1,3,1)',1,[]);
    else
        xOff = 0;
        yOff = 0;
    end
    
    % Loop through all masks:
    floeMasks = cell(M,L);
    parfor n = 1:numel(floeXYT)
        % Get normalized vertices for floe polygon (+ periodic multiples)
        vnX = [ (floeXYT{n}(1,:).'-Xf)./dx + 1; NaN ] + xOff ;
        vnY = [ (floeXYT{n}(2,:).'-Yf)./dy + 1; NaN ] + yOff ;
        % Only retain floe multiples if they contain a vertex point in domain
        cs = any( vnX>=1 & vnX<=nx & vnY>=1 & vnY<=ny ,1);
        vnX = vnX(:,cs); vnY = vnY(:,cs);
        % Create combined [X,Y] vertex matrix and binanarize
        vertNormN = [vnX(:),vnY(:)];
        mask = getMask( vertNormN,[nx,ny] );
%         mask = mpoly2mask(vertNormN,[ny,nx]); 
        floeMasks{n} = sparse( mask );
    end
    
    
    % Add to structure
    floeOut.floes.floeMasks = floeMasks;

end


function mask = getMask( vert,nxy )

    % Make empty mask grid
    mask = false( nxy(2),nxy(1) );

    % Extract subset of domain and create mask
    % (faster than calling poly2mask on whole domain [?])
    extMin = max( [floor( min(vert));[0,0]] );
    extMax = min( [ceil(max(vert));nxy] );
    nc = extMax-extMin;
    vertN = (vert-extMin);
    maskN = mpoly2mask(vertN,fliplr(nc));

    % Write back to main grid
    xInd = extMin(1)+(1:nc(1));
    yInd = extMin(2)+(1:nc(2));
    mask(yInd,xInd) = maskN;
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
