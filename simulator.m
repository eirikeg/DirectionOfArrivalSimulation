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
fc = 5.8e9;                                                               %Carrier frequency of signal [Hz]
c = physconst('LightSpeed'); 
lambda = c/fc;                                                            %wavelength [m]

%% Create antenna element on UAV, platform and radiator
antenna = phased.IsotropicAntennaElement;
antennaTX = antenna;

%%%
%TRAJECTORY TO BE CHANGED... Figure out how to generate a feasible
%trajectory.
%%%
%trajectory specification for UAV - from example:  https://se.mathworks.com/help/phased/ref/phased.platform-system-object.html
%piecewise cubic interpolation on the waypoints to derive the position and velocity at each time step

x0 = 5;
y0 = 0;
z0 = 0;
vx = 5;
vy = 10;
vz = 0;
ax = 1;
ay = -1;

t = [0:10];
x = x0 + vx*t + ax/2*t.^2;
y = y0 + vy*t + ay/2*t.^2;
z = z0 + t;
waypoints = [t.' x.' y.' z.'];


antennaOrientation = roty(0);
UAV = phased.Platform('MotionModel','Custom','CustomTrajectory',waypoints,...
                                'InitialOrientationAxes', antennaOrientation,...
                                'OrientationAxesOutputPort',true);
                            
radiator = phased.Radiator('Sensor',antennaTX,'PropagationSpeed',c,...     %Narrowband signal radiator: COnverts signal into radiated wave fields
    'WeightsInputPort',false,'OperatingFrequency',fc);


%% Create stationary phased array antenna ground station
%default propagation speed = speed of light

initPosPARS = [0,0,0]';                                                    %Initial position [m]
initVelPARS = [0,0,0]';                                                    %Inital velocity [m]

arrayOrientation = rotz(0);
arrayRX = phased.URA('Element', antenna, 'Size', [10,10], 'ElementSpacing',...
                    [lambda/2, lambda/2]);                                 %spacing is chosen to avoid spacial aliasing
                
phasedArray = phased.Platform('InitialPosition', initPosPARS, 'Velocity', initVelPARS,...
                                'InitialOrientationAxes', arrayOrientation,...
                                'OrientationAxesOutputPort',true);

%% Transmitter, freeSpace channel, signal

transmitter = phased.Transmitter('PeakPower',1000.0,'Gain',40);            % Output power and antenna gain - USED THE ONE FROM EXAMPLE
channel = phased.FreeSpace('OperatingFrequency', fc);                      % Used for signal propagation from one point to another


%%%
%UNSURE ABOUT HOW TO GENERATE SIGNAL? IS A SIMPLE SINEWAVE CORRECT?
%%%
tau = 100e-6;                                                              %Pulsewidth
prf = 5000;                                                                %Pulse repetition frequency
fs = 1e6;                                                                  %sample rate [Hz]
waveform = phased.LinearFMWaveform('PulseWidth',tau,...                    %signal waveform
    'OutputFormat','Pulses','NumPulses',1,'PRF',prf,'SampleRate',fs);


%Complex sine wave
 dt = 1e-11;
 t = 0:dt:(100*dt);
 w = 5.8e9*2*pi;
 xx = exp(j*w*t);

%% DOA estimator - MUSIC
% 
% estimator = phased.MUSICEstimator2D('SensorArray', arrayRX, 'OperatingFrequency',fc,...
%     'NumSignalsSource','Property',...
%     'DOAOutputPort',true,'NumSignals',1,...
%     'AzimuthScanAngles',-90:1:90,...
%     'ElevationScanAngles',-40:1:40);

%% Simulation loop
stepsize = 0.25;
N = 40;                                                                    %number of steps

azimuth = zeros(1,N);                                                      %Allocate memory for storing data
elevation = zeros(1,N);
% azimuthDOA = zeros(1,N);
% elevationDOA = zeros(1,N);
uavPosition = zeros(3,N);
arrayPosition = zeros(3,N);
timesteps = zeros(1,N);

for t = 1:N+1
    timesteps(t) = t;
    [uavPos ,uavVel] = UAV(stepsize);
    [phasedArrayPos ,phasedArrayVel] = phasedArray(stepsize);
    
    [range,angles] = rangeangle(uavPos, phasedArrayPos, arrayOrientation); %get true propagation path angles and range with respect to the phased array
    
    azimuth(t) = angles(1);                                                %store true azimuth and elevation angle
    elevation(t) = angles(2);
    uavPosition(:,t) = uavPos';
    arrayPosition(:,t) = phasedArrayPos';
    
    signal = xx';                                                          %signal to be transmitted
    radiatedSignal = radiator(signal, angles);
    porpagatedSignal = channel(radiatedSignal,uavPos,phasedArrayPos,...
                                uavVel,phasedArrayVel);
                            
    collectedSignal = collectPlaneWave(arrayRX, porpagatedSignal,...
                                       angles, fc);
%     Run DOA algorithm
%     [~,doas] = estimator(collectedSignal);                                 %DoA Estimation
%     azimuthDOA(t) = doas(1);                                               %Store doa estimates for comparisons
%     elevationDOA(t) = doas(2);
    
end

%% Plot

figure
plot(timesteps, azimuth,'r');
grid on;
hold on;
% plot(timesteps, azimuthDOA,'ro');
plot(timesteps, elevation, 'b');
% plot(timesteps, elevationDOA,'bo');
% legend('azimuth', 'azimuthDOA', 'elevation', 'elevationDOA');
legend('azimuth','elevation');

figure
plot3(x,y,z,'o'); 
hold on;
grid on;
plot3(uavPosition(1,:), uavPosition(2,:), uavPosition(3,:));
plot3(0,0,0,'ro')

