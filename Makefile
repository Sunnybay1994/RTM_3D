fdtdcc=mpicxx
fdtdcflags=
fdtdsrc=src/FDTD
fdtd : $(fdtdsrc)/main.cpp $(fdtdsrc)/updateEH.cpp $(fdtdsrc)/input.cpp $(fdtdsrc)/output.cpp 
	$(fdtdcc) $(fdtdcflags) -o FDTD_MPI.exe  $(fdtdsrc)/fdtd.h $(fdtdsrc)/main.cpp $(fdtdsrc)/updateEH.cpp $(fdtdsrc)/input.cpp $(fdtdsrc)/output.cpp

pstdcc=ifort
pstdcflags=-qopenmp -mkl
pstdsrc=src/PSTD
pstd : $(pstdsrc)/mkl_dfti.f90 $(pstdsrc)/sgpstd3d.f90
	$(pstdcc) $(pstdcflags) -o PSTD.exe $(pstdsrc)/mkl_dfti.f90 $(pstdsrc)/sgpstd3d.f90

clean :
	rm mkl_dft_type.mod mkl_dfti.mod

all: fdtd pstd clean