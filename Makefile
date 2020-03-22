cc=mpicxx
cflags=
FDTD_MPI_geop : main.cpp updateEH.cpp input.cpp output.cpp 
	$(cc) $(cflags) -o FDTD_MPI_geop main.cpp updateEH.cpp input.cpp output.cpp

