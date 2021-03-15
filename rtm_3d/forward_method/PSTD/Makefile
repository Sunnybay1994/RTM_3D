cc=ifort
cflags=-qopenmp -mkl
PSTD : mkl_dfti.f90 sgpstd3d.f90
	$(cc) $(cflags) -o PSTD.exe mkl_dfti.f90 sgpstd3d.f90