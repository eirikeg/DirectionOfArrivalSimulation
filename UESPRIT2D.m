function [DOA] = UESPRIT2D(x, d, lambda, Mx, My)
%Implementation of 2D ESPRIT
%INPUT:
%   x: collected signal
%   N: Number of samples
%   d: Number of signals
%   Mx: Number of arrays along x-axis, Array has size: Mx X My
%   My: Number of arrays along y-axis, Array has size: Mx X My
%OUTPUT:
%   DOA: Estimate of DOA - [azimuth; elevation]

N = length(x(1,:)); %Number of samples
    
M = length(x(:,1)); %Number of total sensors

% %unitary transform of extended data matrix
PIm = flip(eye(M));     %anti diagonal matrices
PIn = flip(eye(N));
xBar = conj(x);
Z = [x PIm*xBar*PIn];   %extended data matrix - forward backward 

%Theorem 5.1 from introduction to DOA estimation
if mod(M,2) == 0
    Qm = unitaryMtx2n(M/2);
else
    n = floor(M/2);
    Qm = unitaryMtx2n_1(n);
end

Q2n = unitaryMtx2n(N);
GammaX = Qm'*Z*Q2n;     %unitary transform

GammaXGammaX = GammaX*GammaX';
[V,D] = eig(GammaXGammaX);
D = diag(D);
[D, index] = sort(D, 'descend');                %Sort in descending order
V = V(:,index);                                 %Sort corresponding eigenvectors in the same orde
Es = V(:, 1:d);                                 %Compute signal subspace Es


%Selection matrices - Maximum overlap
J1_Mx = [eye(Mx-1) zeros(Mx-1,1)];
J2_Mx = [zeros(Mx-1,1) eye(Mx-1)];
J1_My = [zeros(My-1,1) eye(My-1)];
J2_My = [eye(My-1) zeros(My-1,1)];

Jmu1 = kron(eye(My),J1_Mx);
Jmu2 = kron(eye(My),J2_Mx);
Jv1 = kron(J1_My,eye(Mx));
Jv2 = kron(J2_My,eye(Mx));

sizeJmu1 = size(Jmu1,1);
sizeJv1 = size(Jv1,1);

if(mod(sizeJmu1,2) == 0)
    Qmx= unitaryMtx2n(sizeJmu1/2);
else
    Qmx = unitaryMtx2n_1((sizeJmu1-1)/2);
end

if(mod(sizeJv1,2) == 0)
    Qmy= unitaryMtx2n(sizeJv1/2);
else
    Qmy = unitaryMtx2n_1((sizeJv1-1)/2);
end

Kmu1 = Qmx'*(Jmu1+Jmu2)*Qm;
Kmu2 = Qmx'*(Jmu1-Jmu2)*Qm;
Kv1 =  Qmy'*(Jv1+Jv2)*Qm;
Kv2 = Qmy'*(Jv1-Jv2)*Qm;

upsilonMu = linsolve(Kmu1*Es,Kmu2*Es);              %Solve real valued Invariance Equation
upsilonV = linsolve(Kv1*Es, Kv2*Es);

[~, D] = eig(upsilonMu);

%Automatic pairing of spatial frequencies
sumUpsilon =  upsilonMu + j*upsilonV;

[~, sumOmega] = eig(sumUpsilon);

mu = 2*atan(real(sumUpsilon))*lambda/(2*pi*0.5*lambda);
v = 2*atan(imag(sumUpsilon))*lambda/(2*pi*0.5*lambda);

%Obtain DOA estimate
DOA = zeros(2, length(mu));
for i = 1:length(mu);
    DOA(:,i) = [angle(mu(i)-j*v(i)); asin(norm((mu(i)+j*v(i))))];
end
DOA = DOA * (180/pi);

end
