%AVERAGES   Computes the averages of data along the 3rd dimension.
%    Takes into account NaN values and removes outliers (3*std).

function [m,s,e] = averages(data)
if ndims(data)<3
   m = data;
   s = nan(size(data));
   e = s;
   return
elseif ndims(data)>3
   error('wrong input!');
end
N = size(data,3);
n = ones(size(data));
I = isnan(data);
data(I) = 0;
n(I) = 0;
n = sum(n,3);
m = sum(data,3)./n;
s = data-repmat(m,[1,1,N]);
s(I) = 0;
s = sqrt(sum(s.^2,3)./(n-1));
I = I | abs(data-repmat(m,[1,1,N])) > repmat(3*s,[1,1,N]);
n = ones(size(data));
n(I) = 0;
data(I) = 0;
n = sum(n,3);
m = sum(data,3)./n;
s = data-repmat(m,[1,1,N]);
s(I) = 0;
s = sqrt(sum(s.^2,3)./(n-1));
s((n-1)<=0) = nan;
e = s.*student_t(n-1)./sqrt(n);


% test this function:
%   data = randn(2,2,10);data(1,:,:)=nan;data(1,2,5)=1;
%   [m,s,e] = averages(data)
