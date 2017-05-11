libblas='/usr/lib/libblas.so';
if exist(libblas,'file')~=2
    system('sudo aptitude install libblas-dev');
end

% This will compile the code on linux
mex('CFLAGS=-std=c99 -fPIC','-DDEFINEUNIX','-largeArrayDims','-lmwblas','mtimesx.c'); 