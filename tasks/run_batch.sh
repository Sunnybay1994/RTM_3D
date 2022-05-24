for i in 01 02 03 04 05 06 07 08 09 10 11 12;do
	cd lydd0920_disease0${i}_pstd8_0o
	rm -r Input
	rm -r Output
	rm -r STD
	cd ..
done;
for i in 02 03 04 05 06 07 08 09 10 11 12;do
	cd lydd0920_disease0${i}_pstd8_0o
	rm *
	# cp -r ../lydd0920_disease001_pstd8_0o/RTM0 .
	cd ..
done;
# for i in 01 02 03 04 05 06 07 08 09 10 11 12;do
# 	cd lydd0920_disease0${i}_pstd8_0o/RTM0/Input
# 	ls *src.in_0000
# 	mv *src.in_0000 src.in_0000
# 	cd ../../..
# done;
for i in 01 02 03 04 05 06 07 08 09 10 11 12;do
	cd lydd0920_disease0${i}_pstd8_0o/log
	pwd
	sbatch script_0000.sh
	cd ../..
done;