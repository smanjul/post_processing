function [u,v,w] = readfld(filenam,n1,n2,n3)

% This functions returns velocities on a cylindrical grid, but along
% cartesian coordinates
% TBD: transform velocities to cylindrical polar form
% filenam: File from which velocity field is read
% n1 x n2 x n3: azimuthal x radial x axial resolution
% azimuthal is equi-spaced
% radial and axial can be uniform or non-uniform
% axial is the vertical coordinate

tot_num = n1*n2*n3;
fid = fopen(filenam);

tline = fgetl(fid) % Header line 1
tline = fgetl(fid) % Header line 2

x_coord = textscan(fid,'%f',tot_num);
y_coord = textscan(fid,'%f',tot_num);
z_coord = textscan(fid,'%f',tot_num);

clear x_coord y_coord z_coord tline

u_cell = textscan(fid,'%f',tot_num);
u = u_cell{1,1}; % velocity along x direction
clear u_cell

v_cell = textscan(fid,'%f',tot_num);
v = v_cell{1,1}; % velocity along y direction
clear v_cell

w_cell = textscan(fid,'%f',tot_num);
w = w_cell{1,1}; % velocity along z direction
clear w_cell