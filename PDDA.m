function [DOA] = PDDA(arrayPos, x, d, lambda, searchStep)
%Implementation of PDDA algorithm
%INPUT:
%   arrayPos: Positions of antenna elements
%   x:        Recieved signal at each antenna
%   searchStep:    Steplength in search through pseudospectrum - default value of 1;
%   d:        Number of signals
%   lambda:   Wavelength of signals
%OUTPUT:
%   DOA: [azimuth, elevation] - Direction of arrival estimate

    if nargin < 5
        searchStep = 1; 
    end

    epsilon = 0.01;                     %Scalar to avoid singularities

    r = arrayPos;                                           %Position of antenna elements
    K = @(azi, el) 2*pi*(1/lambda)*[sind(azi)*cosd(el);...  %Wave vector
               sind(azi)*sind(el); cosd(azi)];
    a = @(r,k) exp(-j*r'*k);                                %Steering vector    

    N = length(x(1,:));                                     %Number of samples

    h = x(1,:);                                             %Two sub matrices
    H = x(2:end, :);
    
    p = (h*H')/(h*h');                                      %Propagator vector: Cross correlation of time-series from each antenna element
    e = [1 p]';                                             %Add unit to the first row representing correlation with itself
    
    
    thetaSearch = 0:searchStep:180;
    phiSearch = 0:searchStep:90;
    P = zeros(length(phiSearch),length(thetaSearch));
    peakIdx = [1000;1000];
    w = -1e11;                                              %Max value in pseudospectrum
    for theta = 1:length(thetaSearch)
        for phi = 1:length(phiSearch)
            k = K(phi,theta);
            P(phi,theta) = abs(a(r,k)'*e)^2;                %Pseudospectrum
            if P(phi,theta) > w
                w = (P(phi,theta));
                peakIdx = [phi,theta];
            end
        end
    end
    
    Ps = w-P;                                               %Subtract global maximum from other points - makes us obtain nulls in DOA of incomming signals
    Ppdda = (Ps+epsilon).^(-1);                             %Pseudospectrum of PDDA
        
    DOA = peakIdx;
    
    w = max(Ppdda,[],'all');
    
    
%     mesh(thetaSearch,phiSearch,abs(Ppdda)/abs(w));        %Plot pseudospectrum
%     grid on;
    
    
end