function Y = bi(X, varargin)
% A = bi(B, i1, i2, i3, ..., iN)
%
% Broadcast indexing. 
%
% X is an N-dimensional array
% i1 ... iN are arrays giving the indices for each dimension.
% If these are all M-dimensional arrays, the output Y is an M-dimensional
% array, each element of which is given by A at the corresponding members of the
% index arrays
%
% if any of the index arrays have size singleton dimensions, these are
% expanded. 
%
% E.g. if X=[1 2; 3 4; 5 6]
% bi(X, [1 3], [1 2]) = [1 6];
% bi(X, [1 3], [1;2]) = [1 5; 2 6];

ZeroIndexArray = 0; % make array of full size 
for i=2:nargin
    ZeroIndexArray = ZeroIndexArray .* varargin{i-1}; % auto-broadcasting!!!
end

inds = cell(nargin-1,1);
for i=2:nargin
    inds{i-1} = ZeroIndexArray + varargin{i-1}; % auto-broadcasting!!!
end

Y = X(sub2ind(size(X), inds{:}));