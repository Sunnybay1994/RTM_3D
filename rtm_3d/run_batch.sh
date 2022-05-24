for i in 01 02 03 04 05 06 07 08 09 10 11 12;do
	python model_em.py --max_cpu 256 --server freeosc --pstd -m z --np 8 --steps zc --model make_model/lydd0920_disease0$i.mat --no_gen_model -n
done;