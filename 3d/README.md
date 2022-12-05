
    First modify fft_theta_3d_mode.m for your file, the same way as you would have done for the 2d case, e.g., number of points, file names, mode number etc.
    Run the file, and it should generate 99 files in the folder.
    Move all these files to a new folder to keep your current folder tidy, e., by running mv 3dmode*.dat mode5/. Here mode5 is the directory that I have created and moved all the related files there.
    Concatenate all these files to one single file by: cat mode5/3dmode*.dat >> mode5/mode5series.dat. This will append all the planes to a single file.
    Run the file mode_construction.f90, which will generate the Tecplot file: gfortran -o 3d_mode mode_construction.f90.
