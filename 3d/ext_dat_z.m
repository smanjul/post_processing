function [xy_pl] = ext_dat_z(A,n,n1,n2)

% A : Variable from which z plane data is extracted
% n: Index of the z plane from which data is to be extracted
% n <= n3

cnt1 = (n - 1)*n2*n1; % Counter for z plane
% cnt1 gives the number of elements starting from bottom plate till the
% previous plane.
for ky = 1:n2
    cnt_st = cnt1 + (ky-1)*n1 + 1; % Starting index of the ky plane within z plane
    cnt_end = cnt1 + (ky-1)*n1 + n1;
    work = A(cnt_st:cnt_end);
%     whos work
    work = work(1:n1-1); % truncate for FFT
    xy_pl(ky,:) = work;    
end

xy_pl = xy_pl';