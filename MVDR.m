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
K = @(azi, el) 2*pi*(1/lambda)*[sind(azi)*cosd(el);...  %Wave vector
               sind(azi)*sind(el); cosd(azi)];
a = @(r,k) exp(-j*r'*k);                                %Steering vector    

N = length(x(1,:));                                     %Number of samples

Rxx = (1/N)*x*x';                                       %Estimated covariance matrix: (Mx+My x Mx+My)

RxxInv = (Rxx)^(-1);                                    %Compute inverse of Rxx

%Search through MVDR pseudospectrum
phiSearch = -90:searchStep:90;
thetaSearch = -90:searchStep:90;
peak_val = 0;
peak_idx = [1000;1000];
p = zeros(181,181);
for phi = 1:length(phiSearch)
    for theta = 1:length(thetaSearch)
        k = K(thetaSearch(theta),phiSearch(phi));
        aa = a(r,k);
        p(theta,phi) = 1/(aa'*RxxInv*aa);
        if abs(p(theta,phi)) > peak_val
            peak_val = abs(p(theta,phi));
            peak_idx = [theta;phi];
        end
    end
end

DOA = peak_idx-90;

% mesh(phiSearch,thetaSearch,abs(p)/abs(peak_val))    %Plot spectrum
% % colorbar;
% % close all;

end
