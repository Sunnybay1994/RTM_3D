destpath=bin

fdtdcc=mpicxx
fdtdcflags=-O3 -march=native
fdtdsrc=rtm_3d/forward_method/FDTD
fdtd : $(fdtdsrc)/main.cpp $(fdtdsrc)/updateEH.cpp $(fdtdsrc)/input.cpp $(fdtdsrc)/output.cpp 
	$(fdtdcc) $(fdtdcflags) -o $(destpath)FDTD_MPI.exe  $(fdtdsrc)/fdtd.h $(fdtdsrc)/main.cpp $(fdtdsrc)/updateEH.cpp $(fdtdsrc)/input.cpp $(fdtdsrc)/output.cpp

pstdcc=ifort
pstdcflags=-qopenmp -mkl
pstdsrc=rtm_3d/forward_method/PSTD
pstd : $(pstdsrc)/mkl_dfti.f90 $(pstdsrc)/sgpstd3d.f90
	$(pstdcc) $(pstdcflags) -o $(destpath)/PSTD.exe $(pstdsrc)/mkl_dfti.f90 $(pstdsrc)/sgpstd3d.f90

clean :
	rm mkl_dft_type.mod mkl_dfti.mod

all: fdtd pstd clean