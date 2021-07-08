function [x,r,r_norm] = omp_specify_atoms(A,b,atom_order)

% Orthogonal Matching Pursuit (OMP)
% 
% Input:
%   A: dictionary (matrix)
%   b: signal 
%   atom_order: atoms selected that will be added.
% Output:
%   x: coeff vector for sparse representation
%   r: residual
%   r_norm: norm of residual

[N,K] = size(A); % N:dim of signal, K:#atoms in dictionary
if (N ~= size(b))
    error('Dimension not matched');
end

S = length(atom_order);     % sparsity level
x = zeros(K,1);      % coefficient (output)
r = b;               % residual of b
omega = zeros(S,1);  % selected support
A_omega = [];        % corresponding columns of A
cnt = 0;
while (cnt < S)  % choose S atoms
    cnt = cnt+1;
    x_tmp = zeros(K,1);
    inds = setdiff([1:K],omega); % iterate all columns except for the chosen ones
    for i = inds
        x_tmp(i) = A(:,i)' * r / norm(A(:,i)); % sol of min ||a'x-b||
    end
    ichosen = atom_order(cnt);
    omega(cnt) = ichosen;
    A_omega = [A_omega A(:,ichosen)];
    x_ls = A_omega \ b;  % Aomega * x_ls = b
    r = b - A_omega * x_ls; % update r
end

for i = 1:S
    x(omega(i)) = x_ls(i); %x_sparse(i).value;
end
r_norm = norm(r);
end
