for (( i = 2; i <=8 ; i++));do
	for (( j = 1; j <=3; j++));do
		cd 'v5${i}z$j_pstd8/log'
		pwd
		bash script_0000.sh
		cd ../..
	done;
done;