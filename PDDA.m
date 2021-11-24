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
    K = @(azi, el) 2*pi*(1/lambda)*[sind(el)*cosd(azi);...  %Wave vector
               sind(azi)*sind(el); cosd(el)];
    a = @(r,k) exp(-j*r'*k);                                %Steering vector  
    
%     a = @(az, el) exp( -1i*2*pi*( Rzyx(0,0,deg2rad(az))*Rzyx(0,deg2rad(el),0)*[0 0 1]' )'*r/lambda )';



    N = length(x(1,:));                                     %Number of samples

    h = x(1,:);                                             %Two sub matrices
    H = x(2:end, :);
    
    p = (h*H')/(h*h');                                      %Propagator vector: Cross correlation of time-series from each antenna element
    e = [1 p]';                                             %Add unit to the first row representing correlation with itself
    
    
    phiSearch = 1:searchStep:360;
    thetaSearch = 1:searchStep:90;
    phiLength = length(phiSearch);
    thetaLength = length(thetaSearch);
    peak_val = 0;
    peak_idx = [inf;inf];
    P = zeros(phiLength,thetaLength);
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
        
    DOA = peakIdx*searchStep;
    
    w = max(Ppdda,[],'all');
    
    mesh(thetaSearch,phiSearch,abs(Ppdda)/abs(w));        %Plot pseudospectrum
%     grid on;
    
    
end