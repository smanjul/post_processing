% This code reads the tecplot file for the VE flow
% and does FFT in the theta direction for a specified
% z-plane.
% The input data file (the tecplot file) is already
% transformed to the Cartesian coordinates.
%
% The wavenumber to be extracted is specified
% in the variable mod_no. and is written in a file
% This output file can be read in an accompanied Fortran
% code to write the data in a tecplot file for visualization.
%-------------------------------------------------------------

filenam = 'medie_tec_5000a2_5.data';

n1 = 385; n2 =197; n3 = 385; 
tot_num = n1*n2*n3;

% Read velocity field
[u,v,w] = readfld(filenam,n1,n2,n3);

pln_no = 10; % z-plane number
% Extract data for a particular plane
[w_eg] = ext_dat_z(w,pln_no,n1,n2);

%% fft for the data
w_amp = fft(w_eg);
w_amp = w_amp/(n1-1);
% Current order w_amp are ordered with k -> 0 -> (n1 - 1)/2 -> -(n1-1)/2 + 1-> -1
% w_amp_shft = fftshift(w_amp)
% w_amp_shft order -(n1-1)/2 + 1 ->  0 -> (n1 - 1)/2

mod_no = 10;
thet_ran = linspace(0,2*pi,n1);
thet_ph = exp(1i*mod_no*thet_ran);
rnd = conj(thet_ph');

% k = mod_no;
k_pos = mod_no + 1;
k_neg = n1 - mod_no;
amp_pk = w_amp(k_pos,:); % amplitude corresponding to k. has dimensions of r
amp_nk = w_amp(k_neg,:); % amplitude corresponding to -k.

mod_phys = real(rnd*amp_pk + conj(rnd)*amp_nk);
dlmwrite('mode10_re5000a2_5.dat', mod_phys, ' ');

r_ran = linspace(0,0.4,n2); % ad hoc for now
[XX,YY] = meshgrid(thet_ran,r_ran);
%whos XX YY
XX = XX'; YY = YY';
%whos XX YY
figure(2);
contourf(XX,YY,mod_phys)

% Polar coordinates to cartesian
%[x,y] = pol2cart(theta_ran,r_ran);
%[x_alt,y_alt] = pol2cart(theta_ran(1:n1-1),r_ran); % Skipping theta = 2pi

