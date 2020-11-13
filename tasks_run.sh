cd tasks
# for((i=16;i>0;i--));
# do
# echo fdtd-$i
# cd pstdtest_800MHz_2.0m_0.5m_fdtd_$i
# mpiexec -np $i ./FDTD_MPI 0 >out.log
# cd ..
# echo pstd-$i
# cd pstdtest_800MHz_2.0m_0.5m_pstd_$i
# ./PSTD 0000 >out.log
# cd ..
# done


tail -3 pstdtest_800MHz_2.0m_0.5m_fdtd_1/out.log > time_fdtd.txt
tail -1 pstdtest_800MHz_2.0m_0.5m_pstd_1/out.log > time_pstd.txt
for((i=2;i<17;i++));
do
tail -3 pstdtest_800MHz_2.0m_0.5m_fdtd_$i/out.log >> time_fdtd.txt
tail -1 pstdtest_800MHz_2.0m_0.5m_pstd_$i/out.log >> time_pstd.txt
done