#!/bin/bash

# for((i=1;i<17;i++));
# do
# echo gen-tasks-$i
# python model_em.py -f 800 --dx_src 2 --dx_rec 0.5 --np $i --half_span 0
# done
# python model_em.py -f 800 --dx_src 2 --dx_rec 0.5 --pstd --half_span 0

cd tasks
for((i=16;i>0;i--));
do
echo run-pstd-$i
cd pstdtest_800MHz_2.0m_0.5m_pstd
./PSTD $i 0 > $i.out
cd ..
echo run-fdtd-$i
cd pstdtest_800MHz_2.0m_0.5m_fdtd_$i
mpiexec -np $i ../../FDTD_MPI 0 > out.log
cd ..
done

tail -3 pstdtest_800MHz_2.0m_0.5m_fdtd_1/out.log > time_fdtd.txt
tail -1 pstdtest_800MHz_2.0m_0.5m_pstd/1.out > time_pstd.txt
for((i=2;i<17;i++));
do
tail -3 pstdtest_800MHz_2.0m_0.5m_fdtd_$i/out.log >> time_fdtd.txt
tail -1 pstdtest_800MHz_2.0m_0.5m_pstd/$i.out >> time_pstd.txt
done