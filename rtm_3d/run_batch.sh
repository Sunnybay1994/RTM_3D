# python model_em.py --model make_model/IF.mat --np 8 --server freeosc --pstd
# for v in -5 -2 -1 1 2 5 ;do
# 		python model_em.py --model make_model/IFv${v}.mat --np 8 --server freeosc --pstd
# done;
# for v in -5 5;do
# 	for z in 1 3;do
# 			python model_em.py --model make_model/IFv${v}z${z}.mat --np 8 --server freeosc --pstd
# 	done;
# done;
for snr in 80 90 100 110 120 ;do
		python model_em.py --model make_model/IF.mat --np 8 --server freeosc --pstd --steps pbic --suffix snr$snr --no_gen_model -n
done;