function [Eig, EigVal] = get_eig(Matrix)
%% returns eigenvectors and eigenvalues with eigenvalues descending.
%Eig(:,n) is the nth eigenvalue of Matrix
[Eig,EigVal] = eig(Matrix,'vector');
[EigVal, EigInd] = sort(EigVal,'descend');
Eig = Eig(:, EigInd);
end

