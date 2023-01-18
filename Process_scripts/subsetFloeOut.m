function floeOut = subsetFloeOut( floeOut, LIA )
% SUBSETFLOEOUT subsets all fields in floeOut (in the time dimension) based
% on provided indices
%   
%   floeOut = subsetFloeOut( floeOut, LIA )
%
%   S.D.Brenner 2022

L = length(floeOut.state.time);

flds1 = fields(floeOut);
for n = 2:length(flds1)
    if isstruct( floeOut.(flds1{n}) )
        flds2 = fields( floeOut.(flds1{n}) );
        for m = 1:length(flds2)
            if size( floeOut.(flds1{n}).(flds2{m}),2 )==L
                floeOut.(flds1{n}).(flds2{m}) = floeOut.(flds1{n}).(flds2{m})(:,LIA);
            end
        end
    else
        continue;
    end
end


end