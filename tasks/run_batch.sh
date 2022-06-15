# for i in 01 02 03 04 05 06 07 08 09 10 11 12;do
# 	cd lydd0920_disease0${i}_pstd8_0o
# 	rm -r Input
# 	rm -r Output
# 	rm -r STD
# 	cd ..
# done;
# for i in 02 03 04 05 06 07 08 09 10 11 12;do
# 	cd lydd0920_disease0${i}_pstd8_0o
# 	rm *
# 	cp -r ../lydd0920_disease001_pstd8_0o/RTM0 .
# 	cd ..
# done;
# for i in 01 02 03 04 05 06 07 08 09 10 11 12;do
# 	cd lydd0920_disease0${i}_pstd8_0o
# 	mv ../lydd0920_disease0${i}_src.in_0000 RTM0/Input/src.in_0000
# 	cd ..
# done;
# for i in 01 02 03 04 05 06 07 08 09 10 11 12;do
# 	cd lydd0920_disease0${i}_pstd8_0o/log
# 	pwd
# 	sbatch script_0000.sh
# 	cd ../..
# done;
# for i in 30 45 60;do
# 	remotedir=/home/mabw/win/run/RTM_3D/tasks/l3f/
# 	fn=l3f${i}15_pstd8
# 	fn0=l3f${i}15z_pstd8_0o
# 	cmd="mkdir $remotedir$fn;mkdir $remotedir$fn0;mkdir $remotedir$fn0/RTM0;mkdir $remotedir$fn0/RTM0/Output;"
# 	ssh mypc $cmd
# 	scp $fn/* mypc:${remotedir}$fn
# 	scp $fn0/* mypc:${remotedir}$fn0
# 	scp -r $fn/Result mypc:${remotedir}$fn
# 	scp $fn0/RTM0/Output/*00_06*.bin mypc:${remotedir}$fn0/RTM0/Output
# done;