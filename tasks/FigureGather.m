%% load files and parameters
homedir = 'MOtest3_pstd6_std';
% load(fullfile(homedir,'pdata_agc.mat'))
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
workdir = fullfile(homedir,'.');
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
xx = unique(recinfo(1,:));
y = unique(recinfo(2,:));
tt = (1:nt)*dt_ori/1e-9;
if ~exist(gatherf,'file')
    gfid = fopen(fullfile(outputdir,'merge_gather_0000.bin'),'r');
    gatherData = fread(gfid,[nt,nrec],'float');
    fclose(gfid);
    gatherData_agc = agc(gatherData);
% for i = 1:length(y)
%     ind{i} = ((1:length(xx))-1)*length(y) + i;
%     gather{i} = gatherData(:,ind{i});
%     gather_agc{i} = gatherData_agc(:,ind{i});
% %     t_mean = mean(gather{i},'all');
% %     t_std = std(gather{i},0,'all');
% %     logic_tr_b = abs(gather{i}) > 6 * t_std;
% %     gather{i}(logic_tr_b) = 0;
%     figure(i);imagesc(gather{i})
%     figure(length(y)+i);imagesc(gather_agc{i})
%     pause(0.1)
% end
%     save(gatherf,'pre_rtm_gathers','ntr','nt','gather_infos','-v7.3')
else
    load(gatherf)
end

figure(1);clf;
set(gcf,'Unit','centimeters')
set(gcf,'Position',[0,0,29.7,21])
subplot(2,1,1)
imagesc(gatherData)
subplot(2,1,2)
imagesc(gatherData_agc)
%%
snr = 100;
gatherData_temp = gatherData;
t_std = std(gatherData_temp,0,'all');
logic_tr_b = abs(gatherData_temp) > t_std;
gatherData_temp(logic_tr_b) = 0;
figure(10);imagesc(gatherData_temp)
gatherData_noise = awgn(gatherData_temp,snr) - gatherData_temp + gatherData;
gatherData_noise_agc = agc(gatherData_noise);
figure(2);clf;
set(gcf,'Unit','centimeters')
set(gcf,'Position',[0,0,29.7,21])
subplot(2,1,1)
imagesc(gatherData_noise)
subplot(2,1,2)
imagesc(gatherData_noise_agc)

gfid = fopen(fullfile(outputdir,sprintf('merge_gather_0000_snr%d.bin',snr)),'w');
fwrite(gfid,gatherData_noise_agc,'float');
fclose(gfid);

%%
gfid = fopen(fullfile(outputdir,'merge_gather_0000_AGC.bin'),'w');
fwrite(gfid,gatherData_agc,'float');
fclose(gfid);
gatherData0 = zeros(size(gatherData));
gfid = fopen(fullfile(outputdir,'merge_gather_0000_ZERO.bin'),'w');
fwrite(gfid,gatherData0,'float');
fclose(gfid);
