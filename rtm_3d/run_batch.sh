for (( i = 2; i <=8 ; i++));do
	for (( j = 1; j <=3; j++));do
		python model_em.py --model make_model/v5${i}z$j.mat --np 8 --max_cpu 64 --server freeosc --steps gfpbic
	done;
done;