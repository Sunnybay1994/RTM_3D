cc=ifort
cflags=-qopenmp -mkl -O3 -march=skylake-avx512
PSTD : mkl_dfti.f90 sgpstd3d.f90
	$(cc) $(cflags) -o PSTD.exe mkl_dfti.f90 sgpstd3d.f90