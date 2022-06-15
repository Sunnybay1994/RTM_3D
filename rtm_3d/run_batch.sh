for i in 30 45 60;do
	python model_em.py --max_cpu 64 --server freeosc --pstd --np 8 --model make_model/l3f${i}15.mat --no_gen_model -n
	python model_em.py --max_cpu 128 --server freeosc --pstd -m z --np 8 --model make_model/l3f${i}15z.mat --no_gen_model -n
done;