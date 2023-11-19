destpath=bin/
march=native
# change it according to your CPU may bring a better performance, such as:
# freeosc: skylake-avx512
# 6800k: broadwell

all: march fdtd pstd clean

march:
	@echo march=${march}

fdtdcc=mpicxx
fdtdcflags=-O3 -march=${march}
fdtdsrc=rtm_3d/forward_method/FDTD
fdtd : $(fdtdsrc)/main.cpp $(fdtdsrc)/updateEH.cpp $(fdtdsrc)/input.cpp $(fdtdsrc)/output.cpp 
	$(fdtdcc) $(fdtdcflags) -o $(destpath)FDTD_MPI.exe  $(fdtdsrc)/fdtd.h $(fdtdsrc)/main.cpp $(fdtdsrc)/updateEH.cpp $(fdtdsrc)/input.cpp $(fdtdsrc)/output.cpp

pstdcc=ifort
pstdcflags=-qopenmp -qmkl -O3 -march=${march}
pstdsrc=rtm_3d/forward_method/PSTD
pstd : $(pstdsrc)/mkl_dfti.f90 $(pstdsrc)/sgpstd3d.f90
	$(pstdcc) $(pstdcflags) -o $(destpath)PSTD.exe $(pstdsrc)/mkl_dfti.f90 $(pstdsrc)/sgpstd3d.f90

clean :
	rm mkl_dft_type.mod mkl_dfti.mod