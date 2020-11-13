for((i=1;i<17;i++));
do
echo $i
python model_em.py -f 800 --dx_src 2 --dx_rec 0.5 --np $i --half_span 0
python model_em.py -f 800 --dx_src 2 --dx_rec 0.5 --nthrd $i --pstd --half_span 0
done