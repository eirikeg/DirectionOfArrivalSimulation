function [DOA] = ESPRIT(arrayPos, x, d, lambda)
%Implementation of 2D ESPRIT
%arrayPos: Position of array
%x: collected signal
%N: Number of samples
%d: Number of signals
c = physconst('LightSpeed');
f = 1e9;
r = arrayPos;
K = @(azi, el) 2*pi*(1/lambda)*[sind(azi)*cosd(el); sind(azi)*sind(el); cosd(azi)];

M = length(x(:,1));
N = length(x(1,:));

R1 = x(1:M-1,:);
R2 = x(2:M,:);

R = [R1; R2];
Rrr = (1/N)*(R*R');

[V, D] = eig(Rrr);
D = diag(D);                                     %Create vector of the eigenvalues
[D, index] = sort(D, 'descend');          %Sort in descending order
V = V(:,index);                                 %Sort corresponding eigenvectors in the same orde
Vs = V(:, 1:d);

Vs1 = Vs(1:(M-1),:);
Vs2 = Vs(M:end,:);

C = [Vs1';Vs2']*[Vs1 Vs2];

[Vc, Dc] = eig(C);
Dc = diag(Dc);                                     %Create vector of the eigenvalues
[Dc, index] = sort(Dc, 'descend');          %Sort in descending order
Vc = Vc(:,index);                                 %Sort corresponding eigenvectors in the same order

psi = -Vc(1,2)*Vc(2,2)^(-1);
[~, Dpsi] = eig(psi);
Dpsi = diag(Dpsi);


DOA = asin(c/(2*pi*f*(lambda*0.5))*angle(max(Dpsi)))*(180/pi);        

end
