function v = IndexArrayNan(a, i)
% v = IndexArrayNan(a, i);
% 
% given an n-dim array a, and an array i of size (n, d1, ...) array i, return
% all values of a given by the indices in i, returning NaN if out of range
%
% so v will be of size (d1, ...)
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
 

sza = size(a);
szi = size(i);
szv = szi(2:end); % size of output array

% make i into a n by D=d1*d2*... matrix
i2 = i(:,:);

% get number of dimensions in actual array and index array
nda = length(sza);
ndi = szi(1);

% AllInRange is a 1xD matrix saying whether all coordinates are in range
if ndi==nda
    CoordsInRange = (i2>=1 & bsxfun(@le, i2, sza'));
    AllInRange = all(CoordsInRange,1);
elseif ndi>nda
    % if they are not equal, check that all trailing dimensions of index array are 1
    CoordsInRange = (i2(1:nda,:)>=1 & bsxfun(@le, i2(1:nda,:), sza'));
    ExtrasOne = (i2(nda+1:end,:)==1);
    AllInRange = all(CoordsInRange,1) & all(ExtrasOne,1);
else
    error('Index array not enough dimensions');
end
   

% D size row vector

% create linear index from i2. (God i hate matlab for having to do this cell thing)
LinIndCell = num2cell(i2(:,AllInRange),2);

% create output
if length(szv)==1 % now i REALLY hate matlab.
    v = nan(szv,1); % because otherwise it returns a square
else
    v = nan(szv);
end

if iscell(a)
    v = num2cell(v);
end

v(AllInRange) = a(sub2ind(sza, LinIndCell{:}));
return
