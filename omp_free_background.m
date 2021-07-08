function [x,r,r_norm] = omp_free_background(A,b,n_nonzero_coefs,thresh,background_indices)

% Orthogonal Matching Pursuit (OMP)
% https://github.com/seunghwanyoo/omp
% Input:
%   A: dictionary (matrix)
%   b: signal 
%   n_nonzero_coefs: max number of atoms from A (not including background)
%                    that can be selected.
%   thresh: OMP stops when reduction in residue drops below thresh.
%   background_indices: atoms at these indices will be selected first.
% Output:
%   x: coeff vector for sparse representation
%   r: residual
%   r_norm: norm of residual

[N,K] = size(A); % N:dim of signal, K:#atoms in dictionary
if (N ~= size(b))
    error('Dimension not matched');
end

if nargin<3 || isempty(n_nonzero_coefs)
    n_nonzero_coefs=round(K/10);
end
if nargin<4 || isempty(thresh)
    thresh=0;
end
if nargin<5 || isempty(background_indices)
    background_indices=[];
end


S_background = length(background_indices);
S = S_background+n_nonzero_coefs;     % sparsity level
x = zeros(K,1);      % coefficient (output)
r = b;               % residual of b
r_norm = norm(b);    % norm of residual of b
r_norm_last = inf;   % inf so have to have atleast one coefficient.
omega = zeros(S,1);  % selected support
omega_last = omega;
A_omega = [];        % corresponding columns of A
cnt = 0;
while (cnt < S)  % choose S atoms
    cnt = cnt+1;
    if cnt<=S_background
        ichosen = background_indices(cnt);
    else
        x_tmp = zeros(K,1);
        inds = quick_setdiff(1:K,omega(omega~=0)); % iterate all columns except for the chosen ones
        for i = inds
            x_tmp(i) = A(:,i)' * r / norm(A(:,i)); % sol of min ||a'x-b||
        end
        [~,ichosen] = max(abs(x_tmp)); % choose the maximum
    end
    omega(cnt) = ichosen;
    A_omega = [A_omega A(:,ichosen)];
    x_ls = A_omega \ b;  % Aomega * x_ls = b
    r = b - A_omega * x_ls; % update r
    r_norm = norm(r);
    if r_norm_last-r_norm<thresh && cnt>S_background
        %If met residue threshold, go back to last solution 
        %i.e. current solution is too many atoms.
        x_ls = x_ls_last;
        omega = omega_last(1:cnt-1);
        r_norm = r_norm_last;
        S = cnt-1;
        break;
    else
        x_ls_last = x_ls;
        omega_last = omega;
        r_norm_last = r_norm;
    end
end

for i = 1:S
    x(omega(i)) = x_ls(i); %x_sparse(i).value;
end
end

function Z = quick_setdiff(X,Y)
if ~isempty(X)&&~isempty(Y)
    check = false(1, max(max(X), max(Y)));
    check(X) = true;
    check(Y) = false;
    Z = X(check(X));
else
    Z = X;
end
end
