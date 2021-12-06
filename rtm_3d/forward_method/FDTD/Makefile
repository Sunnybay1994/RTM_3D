cc=mpicxx
cflags=-O3 -march=skylake-avx512
FDTD_MPI : main.cpp updateEH.cpp input.cpp output.cpp 
	$(cc) $(cflags) -o FDTD_MPI.exe  fdtd.h main.cpp updateEH.cpp input.cpp output.cpp

