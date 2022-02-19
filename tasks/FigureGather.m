%% load files and parameters
homedir = 'cntc_400MHz_5x5_0_0.1x0.5_pstd_6';
% srcpulse_type = 'ricker';
srcpulse_type = 'ricker';
it0_ori = 0;
if strcmp(srcpulse_type,'ricker')
    if contains(homedir,'100MHz')
        it0_ori = 401;
    elseif contains(homedir,'400MHz')
        it0_ori = 101;
    end
end
% workdir = fullfile(homedir);
workdir = fullfile(homedir,'STD');
inputdir = fullfile(workdir,'Input');
outputdir = fullfile(workdir,'Output');
fxslice = dir(fullfile(outputdir,'slx*.bin'));
fyslice = dir(fullfile(outputdir,'sly*.bin'));
fzslice = dir(fullfile(outputdir,'slz*.bin'));
fwave = dir(fullfile(outputdir,'wvf*.bin'));
outdir = workdir;
slicedir = fullfile(workdir,'slices');
if ~exist(slicedir,'dir')
    mkdir(slicedir)
end

modelfiles = dir(fullfile(homedir,'model.mat'));
modelfn = fullfile(modelfiles.folder,modelfiles.name);
m = load(modelfn);

dx_ori = m.dx;
dy_ori = m.dy;
dz_ori = m.dz;

nx_ori = m.nx;
ny_ori =m.ny;
nz_air_ori = m.nz_air;
nz_ori = m.nz;

dt_ori = m.dt;%ns

slicex = m.slicex;
slicey = m.slicey;
slicez = m.slicez;

slice_outstop = m.outstep_slice;
x_outstep = m.outstep_x_wavefield;
t_outstep = m.outstep_t_wavefield;

srcinfo = m.srcinfo;
recinfo = m.recinfo;
nsrc = size(srcinfo,2);
nrec = size(recinfo,2);

nx = nx_ori / x_outstep;
ny = ny_ori / x_outstep;
nz = nz_ori / x_outstep;
nz_air = nz_air_ori / x_outstep;
dx = dx_ori * x_outstep;
dy = dy_ori * x_outstep;
dz = dz_ori * x_outstep;

dt = dt_ori*t_outstep;
it0 = round(it0_ori / t_outstep);

%% read nt
fid = fopen(fullfile(inputdir,'par.in'),'r');
C = textscan(fid,'%d',4,'headerlines',3,'Delimiter',',');
nt = double(C{1}(4));
fclose(fid);

%% behavier
draw_slices = true;
draw_wyslice = true;
draw_wzslice = true;
draw_3D_view = false;

%% figure gathers
gatherf = fullfile(outputdir,'gather.mat');
xx = 0:0.1:8;
xx0 = [6:dx:9.6,10.3:dx:14.2]-6;
tt = (1:nt)*dt_ori/1e-9;
if ~exist(gatherf,'file')
    gfid = fopen(fullfile(outputdir,'gather_0000_00000.bin'),'r');
    gatherData = fread(gfid,[nt,nrec],'float');
    fclose(gfid);
    ntr = [81,81,81,81,77];
    figure(1);gather{1} = gatherData(:,1:81);
    imagesc(xx,tt,gather{1})
    gather{2} = gatherData(:,81*4+(1:77));
    figure(2);imagesc(xx0,tt,gather{2})
    for i = 3:5
        gather{i} = gatherData(:,81*(i-2)+(1:81));
        figure(i);imagesc(xx,tt,gather{i})
    end
%     save(gatherf,'pre_rtm_gathers','ntr','nt','gather_infos','-v7.3')
else
    load(gatherf)
end
