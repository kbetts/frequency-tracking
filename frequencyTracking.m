%% Project:  Frequency Tracking
%% Kellen Betts   |  kellen.betts@gmail.com
%% Date:  120227  |  Version:  1.0
%% Description: 	 Time-frequency anaylsis of noisy data with moving signature frequency;
%%					 use of FFT, temporal averaging, and Gaussian filtering.

clear all; close all;

%%==========================================================     initializations

load testdata

L = 15; % length of spatial domain
n = 64; % Fourier modes

% spatial domain
x2 = linspace(-L,L,n+1);
x = x2(1:n); y=x; z=x;
[X,Y,Z] = meshgrid(x,y,z);

% Fourier domain
k = (2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[Kx,Ky,Kz] = meshgrid(k,k,k);
[Ksx,Ksy,Ksz] = meshgrid(ks,ks,ks);

% get # of time slices (generalization)
[slices dummy] = size(Undata);

%%================================================================     averaging

% initialize vectors/matrices
Ut = zeros(slices,n^3);
Utave = zeros(n,n,n);

for j=1:slices
	
	% Fourier operation (3D)
	Un(:,:,:) = reshape(Undata(j,:),n,n,n);
	Unt = fftn(Un);
	
	% averaging
	Utave = Utave + Unt;
	
end

Utave = abs(Utave)/slices;

% find center freq
[maxAve, I] = max( reshape(Utave,1,n^3) );
[Iy Ix Iz] = ind2sub(size(Utave),I);
Fx = k(Ix); Fy = k(Iy); Fz = k(Iz);

% plot averaged data (X-Z) to confirm location
subplot(1,2,1), isosurface(Ksx,Ksy,Ksz,( fftshift(Utave)/maxAve ),0.6)
axis([-5 5 -5 5 -5 5]), axis square,
set(gca,'XTick',[-5:1:5])
set(gca,'ZTick',[-5:1:5])
grid on, grid minor, view([0 -L 0]), drawnow
xlabel('Kx'), ylabel('Ky'), zlabel('Kz')

% plot averaged data (Y-Z) to confirm location
subplot(1,2,2), isosurface(Ksx,Ksy,Ksz,( fftshift(Utave)/maxAve ),0.6)
axis([-5 5 -5 5 -5 5]), axis square,
set(gca,'YTick',[-5:1:5])
set(gca,'ZTick',[-5:1:5])
grid on, grid minor, view([-L 0 0]), drawnow
xlabel('Kx'), ylabel('Ky'), zlabel('Kz')

%%================================================================     filtering

% Gaussian filter
tau = 1;
filter = exp(-tau*( (Kx-Fx).^2 + (Ky-Fy).^2 + (Kz-Fz).^2 ));

loc = zeros(slices,3);

% loop through slices
for j=1:slices
	
	% filter slice
	Un(:,:,:) = reshape(Undata(j,:),n,n,n);
	Unt = fftn(Un);
	Uft = Unt.*filter;
	Uf = ifftn(Uft);
	
	% find center coordinates
	[maxVal, I] = max( reshape(Uf,1,n^3) );
	[Iy Ix Iz] = ind2sub(size(Uf),I);
	%loc(j,1) = Ix; loc(j,2) = Iy; loc(j,3) = Iz;
	%loc(j,1) = k(Ix); loc(j,2) = k(Iy); loc(j,3) = k(Iz);
	loc(j,1) = x(Ix); loc(j,2) = y(Iy); loc(j,3) = z(Iz);
	
	figure(2), isosurface(X,Y,Z,( abs(Uf)/abs(maxVal) ),0.8)
	axis([-L L -L L -L L]), grid on, hold on, drawnow
	
end

% plot marble path overlay
figure(2), plot3(loc(:,1),loc(:,2),loc(:,3))
axis([-L L -L L -L L]), grid on, drawnow
xlabel('X'), ylabel('Y'), zlabel('Z')
hold off

% plot marble path alone
figure(3), plot3(loc(:,1),loc(:,2),loc(:,3))
axis([-L L -L L -L L]), grid on, drawnow
xlabel('X'), ylabel('Y'), zlabel('Z')

%%======================================================================     end