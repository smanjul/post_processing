clear; clc; close all;
gam = 1;
probe=load("Air_3400_Ma_01_vel2_combine_dim_less"); % ENTER file to load
T=0.001*gam; % ENTER the sampling time or 1/data rate
sta=1500000; %0 ENTER instant to start
%endu=40000; % ENTER instant upto which signal need to be taken
time=probe(sta:end,1); % ENTER colm which is time
time=gam*time;   % multiplied by aspect ratio to compare with literature

%% Read different columns of the probe data
% First column of the file is time
% Second column is azimuthal velocity
% Third column is radial velocity
% Fourth column is axial velocity
v1=probe(sta:end,2);
v2=probe(sta:end,3);
v3=probe(sta:end,4);
vmag=sqrt(v1.^2+v2.^2+v3.^2);
%vmag=v1;
%vmag=probe(sta:end,5); % ENTER column which is the variable
figure()
plot(time,vmag);
title({'Velocity'},'Interpreter','latex','FontSize',20);
%% Signal
V_mean=mean(vmag);
vmagN=vmag-V_mean;
t = time(1:end);
Fs = 1/T;
figure()
 plot(t,vmagN)
 title({'Mean subtracted Velocity'},'Interpreter','latex','FontSize',20);
 
%% FFT
N = size(vmagN,1);     
NFFT = 2^nextpow2(N);  
P = fft(vmagN, NFFT)/N;
f = Fs/2*linspace(0, 1, NFFT/2+1);
% Power
P1 = 2*abs(P(1:NFFT/2+1)).^2; 
% Amplitude
% P2 = Pmax*2*abs(P(1:NFFT/2+1));  
P2 = 2*abs(P(1:NFFT/2+1));
figure()
%semilogy(f,P2,'k') 
phys_time=1./f;
plot(phys_time,P2,'k')
grid on;
title('Single-sided amplitude spectrum of Velmag(t)','Interpreter','latex','FontSize',20)
xlabel('Frequency (Hz)')
ylabel('|P(f)|')
%+xlim([0 2])
x=vmagN;
%% Autocorrelation
[Rxx,lag]=xcorr(vmag);
RxxN=xcorr(vmagN);
tlag=lag*T;
figure()
plot(tlag,Rxx);
%hold on
%plot(RxxN);

%% Dowsampling
% Interpolate the data to a coarser mesh for a fast calculation of E1 and
% E2 (Cao's method)
%len = length(time);
%xf = linspace(0, 0.001, len);
xf = time(1:end);
yf = x;
figure()
plot(xf,yf)
%%
[xf, index] = unique(xf);            % removing duplicate points
index = sort(index);                 % sorting the index as function 'unique' changes the sorting of index
yf = yf(index);
len = length(xf);
len = floor(0.05*len);                 % only integer values
figure()
plot(xf,yf)
%%
%dx = (time(end)-time(1))/(len-1)
dx = (xf(end)-xf(1))/(len-1);          % to make sure that indices are sorted
%xi = linspace(0, 0.01, len/10);
xi(1)=xf(1);
for i=2:len
  xi(i) = xi(i-1)+dx;
end
%xi=xi';
yi = interp1(xf, yf, xi, 'linear','extrap');
xi=xi';
yi=yi';
%%
xiyi = [xi(:), yi(:)];
dlmwrite('undersampled160.dat', xiyi, 'delimiter', '\t');


 %% Lyapunov
 xdata = yi;%vmagN; %data(:,1);
 dim = 51;
 [~,lag] = phaseSpaceReconstruction(xdata,[],dim);
 eRange = 500;
 lyapunovExponent(xdata,Fs,lag,dim,'ExpansionRange',eRange)
%%
%Lyapunov exponent input is  velocity magnitude  and frequency 
%lyapExp = lyapunovExponent(yi(:),dx);

%%
%figure(1)
%plot(xf, yf, '-b')
%hold on
%plot(xi, yi, 'pg')
%hold off
%% Time lag

%lagf=1000;% ENTER upto what lag in data points in calculate AMI
%lag=[1:lagf];  % 
%v=mai(yi,lag);
%figure()
%plot(lag,v)
%% Embedding dimension
%tao=134; %  ENTER the time delay in data points
tao=750; %  ENTER the time delay in data points
%tao=200; %  ENTER the time delay in data points
%mmax=5; %  ENTER how many dimensions to CHECK
%[E1 E2] = cao_deneme(x,tao,mmax);
%figure()
x1=x(1:end-(2*tao));
x2=x(1+tao:end-tao);
x3=x(1+(2*tao):end);
%sta1=1;       % ENTER start of data to plot
%endu1=100000;     % ENTER max length to plot
%endu1=N-(2*tao)
%plot3(x1(sta1:endu1),x2(sta1:endu1),x3(sta1:endu1))
plot3(x1,x2,x3)


