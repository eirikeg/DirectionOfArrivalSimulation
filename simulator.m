%SIMULATOR FOR PARS UAV SYSTEM
% This script aims to simulate a moving UAV with a single transmitting
% antenna which moves in space. A phased array antenna is located in the origin of the global coordinate system. 
%
% Toolboxes used: Phased Array System Toolbox
% 
% The simulator is a modification of the following examples: 
% https://se.mathworks.com/help/phased/ug/local-and-global-coordinate-systems-example-in-radar.html
%% Clean up
clc; clear;

%% Constants
f = 2.4e9;                                                               %Carrier frequency of signal [Hz]
c = physconst('LightSpeed'); 
lambda = c/f;                                                            %wavelength [m]

%% Create antenna element on UAV, platform and radiator
antenna = phased.IsotropicAntennaElement;
antennaTX = antenna;

%TRAJECTORY for UAV
t = [0:200:800];
xPoints = [400;4000;4000; 1000; 0]';
yPoints = [300; 1000; 4000; 4000; 10]';
zPoints = [0;500; 1000; 500; 0]';

waypoints = [t.' xPoints.' yPoints.' zPoints.'];


antennaOrientation = rotz(0);
UAV = phased.Platform('MotionModel','Custom','CustomTrajectory',waypoints,...
                                'InitialOrientationAxes', antennaOrientation,...
                                'OrientationAxesOutputPort',true);
                            
radiator = phased.Radiator('Sensor',antennaTX,'PropagationSpeed',c,...     %Narrowband signal radiator: COnverts signal into radiated wave fields
    'WeightsInputPort',false,'OperatingFrequency',f);


%% Create stationary phased array antenna ground station
%default propagation speed = speed of light

initPosPARS = [0,0,0]';                                                    %Initial position [m]
initVelPARS = [0,0,0]';                                                    %Inital velocity [m]

Mx = 3;
My = 3;

arrayOrientation = rotz(0);
arrayRX = phased.URA('Element', antenna, 'Size', [Mx,My], 'ElementSpacing',...
                    [lambda/2, lambda/2], 'ArrayNormal', 'z');                                 %spacing is chosen to avoid spacial aliasing
                
phasedArray = phased.Platform('InitialPosition', initPosPARS, 'Velocity', initVelPARS,...
                                'InitialOrientationAxes', arrayOrientation,...
                                'OrientationAxesOutputPort',true);

arrayPos = getElementPosition(arrayRX);
%% Transmitter, freeSpace channel, signal

transmitter = phased.Transmitter('PeakPower',1000.0,'Gain',40);            % Output power and antenna gain - USED THE ONE FROM EXAMPLE
channel = phased.FreeSpace('OperatingFrequency', f);                       % Used for signal propagation from one point to another


%Signal - phasor
fs = 2*f;
t = 0:2e-6:(160-2)*1e-6;
Nsamples = length(t);
s = exp(j*2*pi*f*t);

%% Simulation loop
stepsize = 0.1;
N = 3*(1/stepsize)                                                         %number of steps

azimuth = zeros(1,N);                                                      %Allocate memory for storing data
elevation = zeros(1,N);
azimuthDOA = zeros(1,N);
elevationDOA = zeros(1,N);
uavPosition = zeros(3,N);
arrayPosition = zeros(3,N);
timesteps = zeros(1,N);
uavVelocity = zeros(3,N);
ranges = zeros(1,N);

r = arrayPos;
K = @(azi, el) 2*pi*(1/lambda)*[sind(azi)*cosd(el); sind(azi)*sind(el); cosd(azi)];
a = @(r,k) exp(-j*r'*k);

signal = zeros(1,Nsamples);
for t = 1:N+1
    timesteps(t) = t;
    
    [uavPos ,uavVel, uavOrientation] = UAV(stepsize);
    [phasedArrayPos ,phasedArrayVel] = phasedArray(stepsize);

    %TRASMIT SIGNAL
    %Create UAV signal platform transmitting signal for 160 microseconds
    initPos = uavPos-((160e-6)*uavVel); 
    initVel = uavVel;
    uavSignalPlatform = phased.Platform('MotionModel','Acceleration','InitialPosition', initPos, 'InitialVelocity', initVel,...
                                'InitialOrientationAxes', uavOrientation, 'InitialOrientationAxes', uavOrientation);
    for i=1:80
        [pos, ~] = uavSignalPlatform(i*(1e-6));
        sample = transmitter(s(i));
        sample = channel(sample,pos,phasedArrayPos,...
                                    uavVel,phasedArrayVel);
        signal(i) = sample;
        
    end
    %END OF SIGNAL TRANSMISSION
    
    
    %get true propagation path angles and range with respect to the phased array
    [range,angles] = rangeangle(uavPos, phasedArrayPos, arrayOrientation);
    
    %Store values
    ranges(t) = range;
    azimuth(t) = angles(1);                                                
    elevation(t) = angles(2);
    uavPosition(:,t) = uavPos';
    arrayPosition(:,t) = phasedArrayPos';
    uavVelocity(:,t) = uavVel;
    
    %RECIEVED SIGNAL
    k = K(azimuth(t), elevation(t));    %wavevector

    x = a(r,k)*signal;  
    x = awgn(x, 20);                    %Add gaussian white noise 
    
    %DOA ESTIMATE
    [DOA] = MUSIC(r,x,1,lambda,1);
%     DOA=PDDA(arrayPos, x, 1, lambda, 1);
%     DOA = MVDR(r,x,1,lambda,1);
    
    azimuthDOA(t) = DOA(1);
    elevationDOA(t) = DOA(2);
%                                    
end

%% Plot

timesteps = timesteps*0.1;

figure
plot(timesteps, azimuth,'b', 'LineWidth', 2);
grid on;
hold on;
plot(timesteps, azimuthDOA,'r.');
plot(timesteps, elevation, 'y', 'LineWidth', 2);
plot(timesteps, elevationDOA,'k.');
legend('azimuth', 'azimuthDOA', 'elevation', 'elevationDOA');
% 
figure
plot3(xPoints,yPoints,zPoints,'o'); 
hold on;
grid on;
plot3(uavPosition(1,:), uavPosition(2,:), uavPosition(3,:));
xlabel('x [m]');
ylabel('y [m]');
zlabel('z [m]');
plot3(0,0,0,'ro');
                 

