%FOLDDATA   Folds the embryo measurements, left-right.

function data = folddata(data)
if ndims(data)>3
   error('wrong input!');
end
sz = size(data);
N = sz(1);    % old size
M = N/2+1;    % new size   %%% If this breaks, remember you need an even number of samples around d/v!
sz(1) = 1;
% data(1,:) and data(M,:) are not duplicated.
z = NaN(sz);
data2 = cat(1,z,data(M+1:end,:,:),z);
data = data(1:M,:,:);
cat(3,data,data2);
