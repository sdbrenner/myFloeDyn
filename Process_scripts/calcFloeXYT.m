function floeOut = calcFloeXYT( floeOut )
% FLOEXYT returns the coordinates defining floe outlines
%
%   floeOut = calcFloeXYT( floeOut ), where 'floeOut' was produced by the
%   'readFloeOut' function. The 'floeOut' structure is modified to include
%   a new field:
%       floeOut.floes.floeXYT
%
%   S.D.Brenner, 2022


[M,L] = size( floeOut.state.xa );
xp = floeOut.state.xp;
yp = floeOut.state.yp;
th = floeOut.state.th;
floeShapes = floeOut.floes.floeShapes;

floeXYT = cell(M,L);
% Loop through floes
for l = 1:L
    for m = 1:M
        % create rotation matrix
        R = [ cos(th(m,l)), -sin(th(m,l)) ;
          sin(th(m,l)),  cos(th(m,l))  ];
        % transform floe
        floeXYT{m,l} = ( [floeShapes{m}(1,:);floeShapes{m}(2,:)]'*R' )' + [xp(m,l);yp(m,l)];
    end
end

floeOut.floes.floeXYT = floeXYT;

end