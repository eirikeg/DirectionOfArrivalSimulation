function DOA = MVDR(arrayPos,x, d, lambda, searchStep)
%Implementation of MVDR beamformer
%Input:
%   X: Recieved signal of N snapshots from URA
%   d: Number of signals
%   lambda: Wavelength of signal
%   searchStep: Step size of search through pseudospecturm
%OUTPUT:
%   DOA: Estimated direction of arrival

if nargin < 5
        searchStep = 1; 
end

r = arrayPos;                                           %Position of antenna elements
K = @(azi, el) 2*pi*(1/lambda)*[sind(el)*cosd(azi);...  %Wave vector
               sind(azi)*sind(el); cosd(el)];
a = @(r,k) exp(-j*r'*k);                                %Steering vector    

N = length(x(1,:));                                     %Number of samples

Rxx = (1/N)*x*x';                                       %Estimated covariance matrix: (Mx+My x Mx+My)

RxxInv = (Rxx)^(-1);                                    %Compute inverse of Rxx

%Search through MUSIC pseudospectrum
phiSearch = 1:searchStep:360;
thetaSearch = 1:searchStep:90;
phiLength = length(phiSearch);
thetaLength = length(thetaSearch);
peakVal = 0;
peakIdx = [inf;inf];
P = zeros(phiLength,thetaLength);
w = -1e11;                                              %Max value in pseudospectrum
for theta = 1:length(thetaSearch)
    for phi = 1:length(phiSearch)
        k = K(phi,theta);
        aa = a(r,k);
        P(phi,theta) = 1/(aa'*RxxInv*aa);               %Pseudospectrum
        if P(phi,theta) > peakVal
            peakVal = (P(phi,theta));
            peakIdx = [phi,theta];
        end
    end
end

DOA = (peakIdx*searchStep);

mesh(thetaSearch,phiSearch,abs(P)/abs(peakVal))    %Plot spectrum
xlabel('\phi [deg]');
ylabel('\theta [deg]');
zlabel('Normalized Power')
colorbar;
% % close all;

end
