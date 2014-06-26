%LOOKUPGENOTYPEID  Returns the DePace lab BID genotype ID given a line
%number.
%
%   ID = LOOKUPGENOTYPEID(L)
%
%   L is a line number (must be numeric).
%
% Tara Martin, 2013-01-08


function id = LookupGenotypeID(L)

if isfloat(L)
    %key = [line_number genotype_ID]
    key = [
        204 106
        207 126
        208 127
        209 128
        210 129
        214 137
        240 144
        320 172
        321 173
        322 174
        323 175
        324 176
        325 177
        326 178
        327 179
        328 180
        329 181
        330 182
        356 189
        450 280
        451 281
        452 282
        453 283
        454 284
        455 285
        456 286
        457 287
        458 288
        459 289
        460 290
        ];
    if isempty(find(key(:,1)==L))
        warning('%d not in list of known line numbers.', L)
    else
        id = key(find(key(:,1)==L),2);
    end
else
    warning('%s is not a number.', L)
end

