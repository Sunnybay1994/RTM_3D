#!/bin/bash

min_core=1
max_core=24

if [[ $1 == '1' ]];then
    for((i=$min_core;i<=$max_core;i++));
    do
        echo tasks-$i
        python model_em.py --np $i --server freeosc -n --steps g
        cd tasks/TvsC_800MHz_1.2m_1.2m_fdtd_$i/log
        sh sub.sh
        cd ../../..
    done
elif [[ $1 == '2' ]];then
    ffdtd=tasks/time_fdtd.txt
    fpstd=tasks/time_pstd.txt
    rm $ffdtd
    rm $fpstd
    for((i=$min_core;i<=$max_core;i++));
    do
        tail -3 tasks/TvsC_800MHz_1.2m_1.2m_fdtd_$i/Output/0.out >> $ffdtd
        # tail -1 pstdtest_800MHz_2.0m_0.5m_pstd/$i.out >> $fpstd
    done
else
    echo $1
fi

# for((i=1;i<=$max_core;i++));
# do
# echo gen-tasks-$i
# python model_em.py --np $i --server freeosc -n --steps g > log/$0.log 2>&1 &
# # echo $!
# jobid[$i]=$!
# echo ${jobid[$i]}
# done
# # python model_em.py -f 800 --dx_src 2 --dx_rec 0.5 --pstd --half_span 0

# cd tasks
# for((i=1;i<=$max_core;i++));
# do
# # echo run-pstd-$i
# # cd pstdtest_800MHz_2.0m_0.5m_pstd
# # ./PSTD $i 0 > $i.out
# # cd ..
# wait ${jobid[$i]}
# echo run-fdtd-$i
# cd pstdtest_800MHz_2.0m_0.5m_fdtd_$i/log
# # mpiexec -np $i ../../FDTD_MPI.exe 0 |tee out.log
# sh sub.sh 
# cd ../..
# done

# sh testgrid0.1_300MHz_10.0m_0.2m_pstd/log/sub_0000.sh
# sh testgrid0.1_300MHz_10.0m_0.2m_fdtd_12/log/sub_0000.sh
# sh testgrid0.08_300MHz_10.0m_0.2m_pstd/log/sub_0000.sh
# sh testgrid0.08_300MHz_10.0m_0.2m_fdtd_12/log/sub_0000.sh
# sh testgrid0.05_300MHz_10.0m_0.2m_pstd/log/sub_0000.sh
# sh testgrid0.05_300MHz_10.0m_0.2m_fdtd_12/log/sub_0000.sh
# sh testgrid0.03_300MHz_10.0m_0.2m_pstd/log/sub_0000.sh
# sh testgrid0.03_300MHz_10.0m_0.2m_fdtd_12/log/sub_0000.sh
# sh testgrid0.05_500MHz_10.0m_0.2m_pstd/log/sub_0000.sh
# sh testgrid0.05_500MHz_10.0m_0.2m_fdtd_12/log/sub_0000.sh
# sh testgrid0.05_800MHz_10.0m_0.2m_pstd/log/sub_0000.sh
# sh testgrid0.05_800MHz_10.0m_0.2m_fdtd_12/log/sub_0000.sh

