close all; clear all; clc;
 
theta = [45];
phi = [45];

d = 2;             %nr signals
M = 9;             %nr elements
N = 101;
f = 1e9;
fs = 2*f;
c = physconst('LightSpeed');
lambda = c/f;
t = 0:(1/fs):100*1/fs;
s = exp(j*2*pi*f*t);

array = phased.URA('Size',[3,3],'ElementSpacing',[0.5*lambda 0.5*lambda], 'ArrayNormal', 'z');
r = getElementPosition(array);
K = @(azi, el) 2*pi*(1/lambda)*[sind(azi)*cosd(el); sind(azi)*sind(el); cosd(azi)];

a = @(r,k) exp(-j*r'*k);

k = K(theta, phi);

x = a(r,k)*s;
x = awgn(x, 0);


% [DOA] = PDDA(r, x, 1, lambda);

% DOA = ESPRIT(r, x, 1, lambda);

% DOA = MUSIC(r,x,1,lambda);

DOA = MVDR(r,x,1,lambda);


