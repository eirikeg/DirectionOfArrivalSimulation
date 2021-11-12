% Simulate nordic array and PDDA algorithm

wavelength = 299792458/2480e6;
number_of_sample_sets = 100;
signal_SNR_db = 10;
az_step = 0.5;
el_step = 0.5;
eps = 0.01;

elevation_true = 0;
azimuth_true = 270;

circular_array = false;
if circular_array
    M = 10; % number of elements, placed in a circle
    spacing = 4*0.05; % element spacing. Increasing this above half wavelength can introduce ambiguities
    r = spacing/(2*sin(pi/M)); % radius of array
    ant_pos = r*[cos(0:2*pi/M:2*pi*(1-1/M)); sin(0:2*pi/M:2*pi*(1-1/M)); zeros(1,M)];
    
%     az_ref = 2*pi/M*(11-1);
%     a_ang = @(az, el) exp(-1i*2*pi/wavelength*r*sin(deg2rad(el))*cos(deg2rad(az)-(0:2*pi/M:2*pi*(1-1/M))'));
%     a = @(az, el) exp(-1i*2*pi*(Rzyx(0,0,deg2rad(az))*Rzyx(0,deg2rad(el),0)*[0 0 1]')'*ant_pos/wavelength)';
else
    ant_pos = [-0.00 -0.10 0;
               -0.00 -0.15 0;
               -0.05 -0.15 0;
               -0.10 -0.15 0;
               -0.15 -0.15 0;
               -0.15 -0.10 0;
               -0.15 -0.05 0;
               -0.15 -0.00 0;
               -0.10 -0.00 0;
               -0.05 -0.00 0;
               -0.00 -0.00 0;
               -0.00 -0.05 0]'; % 11 is origin, axes as defined in locator readme
    % ant_pos(:,1:9) = ant_pos(:,1:9)*4; Non-uniform
end
a = @(az, el) exp(-1i*2*pi*(Rzyx(0,0,deg2rad(az))*Rzyx(0,deg2rad(el),0)*[0 0 1]')'*(ant_pos)/wavelength)';


S = complex(1,0)*ones(1,number_of_sample_sets);
% S = complex(randi([0,1],2,number_of_sample_sets)*2-1,0);
A = a(azimuth_true, elevation_true);% exp(-1i*los'*(ant_pos-ant_pos(:,11))/wavelength)';
% Add measurements to create noise
% X = awgn(A*S,signal_SNR_db);
X = A*S;% + normrnd(0,0.4,size(A*S)) + normrnd(0,0.4,size(A*S))*1i;

% Run PDDA
h = X(1,:);
H = X(2:end,:);
p = h*H'/(h*h');
e = [1 p]';
P = @(az,el) abs(a(az,el)'*e)^2;
% Create pseudospectrum
P_spatial = zeros(floor(360/az_step-1), floor(90/el_step-1));
for az_i = 1:floor(360/az_step)
    for el_i = 1:floor(90/el_step)
        az = (az_i-1)*az_step;
        el = (el_i-1)*el_step;
        P_spatial(az_i,el_i) = P(az,el);
    end
end
% Find max value w
w = max(P_spatial,[],'all');
% P_s = @(az,el) w-P(az,el);
P_s = w-P_spatial;
% P_PDDA = @(az,el) 1/(P_s(az,el)+eps);
P_PDDA = (P_s+eps).^-1;


figure(1)
clf
hold on; grid on
surf(0:az_step:(360-az_step), 0:el_step:(90-el_step), P_spatial','EdgeColor','none')
xlabel('Azimuth  (deg)')
ylabel('Elevation (deg)')

figure(2)
clf
hold on; grid on
polarplot3d(flip([P_spatial; P_spatial(1,:)]',2),'polargrid',{9 4},'polardirection','cw');
clf
hold on; grid on
surf(0:az_step:(360-az_step), 0:el_step:(90-el_step), P_s','EdgeColor','none')
xlabel('Azimuth  (deg)')
ylabel('Elevation (deg)')

figure(3)
clf
hold on; grid on
surf(0:az_step:(360-az_step), 0:el_step:(90-el_step), mag2db(P_PDDA'*eps),'EdgeColor','none')
xlabel('Azimuth  (deg)')
ylabel('Elevation (deg)')
scatter3(azimuth_true,elevation_true,0,150,'r+')
scatter3(200,50,0,150,'r+')
     