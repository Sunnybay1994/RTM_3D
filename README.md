# RTM_3D

reverse time migration for ground penetrating radar

## Prerequisites

- **Linux** environment
- **MATLAB** (version M2017 or above), and add the following tools into Matlab path:
  - [CREWES Matlab® Toolbox](https://www.crewes.org/ResearchLinks/FreeSoftware/index.php)
  - [export_fig](https://ww2.mathworks.cn/matlabcentral/fileexchange/23629-export_fig/)
- **C/C++ compiler** and **OpenMP**, install them using the following code:
```shell
sudo apt install build-essential
```
- **[OpenMPI](https://www.open-mpi.org/)/[MPICH](https://www.mpich.org/)**
- **fortran compiler** (better be '[**ifort**](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#fortran)');
- [**Intel® Math Kernel Library (MKL)**](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html);
- **python** (version 3.x or above) and the following packages:
  - **numpy**, **scipy**, **matplotlib**

## Directory Structure

- RTM_3D/
  - rtm_3d/  主要的代码和脚本
    - common/
      -  \_\_init__.py  通用函数，命令行参数分析
      -  normal_moveout.py  NMO代码，没啥用现在
      -  par_RTM.py  读取全局参数
      -  writesource.py  在任务目录下生成发射源信息文件（包含位置和波形），可以将点源展成面源，以减弱旁瓣
    - forward_method/  波场正演源代码
    - image_condition/  
      - corr_RTM_slice_sub.py  对切片应用互相关成像条件
      - corr_RTM_wavefield_sub.py  对整个波场应用互相关成像条件
      - clean.py  清理中间结果
    - make_model/  模型生成目录
    - model_em.py  根据模型生成工作目录及相关数据
    - batchgen.py  生成批处理命令：对每个源生成相应的提交脚本
    - pre_RTM_sub.py  预处理正演数据的接收波形，为计算反传波场作准备
  - tasks/  任务目录，包含对结果作图的脚本
    - figureResult.m  读取多偏移距模式下的RTM结果并将每个源的结果相加作为最终结果，然后作图
    - FigureOutput.m  读取零偏移距模式下的RTM结果并作图


## compile

1. Add the following line into your '*shrc' file（such as '~/.bashrc'） the first time you install ifort or MKL:
```shell
source /opt/intel/oneapi/setvars.sh # sudo install
# or
source ~/intel/oneapi/setvars.sh # user install
```
2. Run the following command in root directory:

```shell
mkdir bin
make all
```

## Usage
 
1. 在make_model目录下用Matlab生成模型参数: *.mat
2. 进入到rtm_3d目录下，运行model_em.py生成工作目录及相应的参数和数据，示例：
```bash
cd rtm_3d
python model_em.py --max_cpu 256 --server freeosc --pstd -m z --np 8 --steps zc --model make_model/[modelname].mat
# 各参数含义可通过‘-h’命令查看
python model_em.py -h
```
3. 进入到task目录下对应的工作目录中的log文件夹下提交运行脚本:
```bash
cd ../tasks/[workdir]/log
bash sub_script.sh
```
4. 在task目录下利用'figure*.m'脚本读取结果并作图。