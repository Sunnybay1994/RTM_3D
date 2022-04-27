# cd IF_pstd8/log
# pwd
# sbatch sub_script.sh
# cd ../..

# for v in -5 -2 -1 1 2 5 ;do
# 	cd IFv${v}_pstd8/log
# 	pwd
# 	sbatch sub_script.sh
# 	cd ../..
# done;

# for v in -5 5;do
# 	for z in 1 3;do
# 		cd IFv${v}z${z}_pstd8/log
# 		pwd
# 		sbatch sub_script.sh
# 		cd ../..
# 	done;
# done;

for snr in 80 90 100 110 120 ;do
	# cp IF/IF_pstd8/Output/merge_gather* IFsnr${snr}_pstd8/Output
	# cp IF/IF_pstd8/STD/Output/merge_gather* IFsnr${snr}_pstd8/STD/Output
	cd IFsnr${snr}_pstd8/log
	pwd
	sbatch sub_script.sh
	cd ../..
done;