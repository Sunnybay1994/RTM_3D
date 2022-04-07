%% load files and parameters
zd = 2;
workdir = sprintf('MOtest2_1
_pstd6');
resultdir = fullfile(workdir,'Result');
f_wavefield_corr = dir(fullfile(resultdir,'result_wavefield_corr*.???'));

fig_result_dir = workdir;%fullfile(workdir,'Result');
%%
modelfiles = dir(fullfile(workdir,'model.mat'));
modelfn = fullfile(modelfiles.folder,modelfiles.name);
m = load(modelfn);
% load(fullfile(workdir,'model_sr.mat'))

dx = m.dx;
dy = m.dy;
dz = m.dz;
dt = m.dt;

nx = m.nx;
ny =m.ny;
nz_air = m.nz_air;
nz = m.nz;

npmlx = m.npmlx;
npmly = m.npmly;
npmlz = m.npmlz;

slices_ix = m.slicex;
slices_iy = m.slicey;
slices_iz = m.slicez;

x = ((1:nx)-npmlx)*dx;
y = ((1:ny)-npmly)*dy;
z = ((1:nz)-nz_air)*dz;

slice_outstop = m.outstep_slice;
x_outstep = m.outstep_x_wavefield;
t_outstep = m.outstep_t_wavefield;

wnx = nx / x_outstep;
wny = ny / x_outstep;
wnz = nz / x_outstep;
wnz_air = nz_air / x_outstep;
wdx = dx * x_outstep;
wdy = dy * x_outstep;
wdz = dz * x_outstep;
wx = x(1:x_outstep:end);
wy = y(1:x_outstep:end);
wz = z(1:x_outstep:end);
wslicex = round(slices_ix / x_outstep);
wslicey = round(slices_iy / x_outstep);
wslicez = round(slices_iz / x_outstep);

%% Draw slices
slice_reload = true;
fig_para = [0.1, 0];
%% model
r1=0.1;
theta = 0:0.1:2*pi;
% xslice
ym1x = y(slices_iy(4)) + r1 * sin(theta);
zm1x = z(slices_iz(2)) + r1 * cos(theta);
% yslice
xm1y = x(slices_ix(1)) + r1 * cos(theta);
zm1y = z(slices_iz(2)) + r1 * sin(theta);
xm2y = x(slices_ix(2)) + r1 * cos(theta);
zm2y = z(slices_iz(2)) + r1 * sin(theta);
xm3y = x(slices_ix(3)) + r1 * cos(theta);
zm3y = z(slices_iz(2)) + r1 * sin(theta);
% zslice
xm1z = x(slices_ix(1)) + r1 * cos(theta);
ym1z = y(slices_iy(4)) + r1 * sin(theta);
xm2z = x(slices_ix(2)) + r1 * cos(theta);
ym2z = y(slices_iy(4)) + r1 * sin(theta);
xm3z = x(slices_ix(3)) + r1 * cos(theta);
ym3z = y(slices_iy(4)) + r1 * sin(theta);

%% Slices
close all
slicey = {[]};
slicey_sum = {[]};
for i =1:length(slices_iy)
    f_slice = dir(fullfile(resultdir,sprintf('result_ycorr_????_%02d.???',i-1)));
    slice_tag = [sprintf("y=%gm",y(slices_iy(i))),"x(m)","depth(m)"];
    [slicey{i},slicey_sum{i}] =  slice_summation(f_slice,x,z,slice_tag,fig_result_dir,slice_reload,fig_para);
    if i == 4       
        % model
        hold on
        plot(xm1y,zm1y,'r--',xm2y,zm2y,'r--',xm3y,zm3y,'r--')
        hold off
        export_fig(fullfile(fig_result_dir,"slice_" + slice_tag(1) + "_with_model.png"),'-transparent')
    end
end
% export_fig(fullfile(fig_result_dir,"slice_" + slice_tag(1) + "_with_model.pdf"),'-transparent')
%%
slicex = {[]};
slicex_sum = {[]};
for i =1:length(slices_ix)
    f_slice = dir(fullfile(resultdir,sprintf('result_xcorr_????_%02d.???',i-1)));
    slice_tag = [sprintf("x=%gm",x(slices_ix(i))),"y(m)","depth(m)"];
    [slicex{i},slicex_sum{i}] =  slice_summation(f_slice,y,z,slice_tag,fig_result_dir,slice_reload,fig_para);
% model
    hold on
    plot(ym1x,zm1x,'r--')
    plot(ym1x,zm1x+0.3,'r--')
    hold off
    export_fig(fullfile(fig_result_dir,"slice_" + slice_tag(1) + "_with_model.png"),'-transparent')
%     export_fig(fullfile(fig_result_dir,"slice_" + slice_tag(1) + "_with_model.pdf"),'-transparent')
end
%%
slicez = {[]};
slicez_sum = {[]};
fig_paraz = [-inf, fig_para(2)];
for i =1:length(slices_iz)
    f_slice = dir(fullfile(resultdir,sprintf('result_zcorr_????_%02d.???',i-1)));
    slice_tag = [sprintf("z=%gm",z(slices_iz(i))),"x(m)","y(m)"];
    [slicez{i},slicez_sum{i}] =  slice_summation(f_slice,x,y,slice_tag,fig_result_dir,slice_reload,fig_paraz);
    if i == 2
        % model
        hold on
        plot(xm1z,ym1z,'r--',xm2z,ym2z,'r--',xm3z,ym3z,'r--')
        hold off
        export_fig(fullfile(fig_result_dir,"slice_" + slice_tag(1) + "_with_model.png"),'-transparent')
    end
end
% export_fig(fullfile(fig_result_dir,"slice_" + slice_tag(1) + "_with_model.pdf"),'-transparent')

return
%% draw slice in one figure
slicey = {[]};
slicey_sum = {[]};
for i =1:length(slices_iy)
    f_slice = dir(fullfile(resultdir,sprintf('result_ycorr_????_%02d.???',i-1)));
    slice_tag = [sprintf("y=%gm",y(slices_iy(i))),"x(m)","depth(m)"];
    [slicey{i},slicey_sum{i}] =  slice_summation(f_slice,x,z,slice_tag,fig_result_dir,slice_reload,fig_para,subplot(2,2,3));
end
hold on
plot(xm1y,zm1y,'r--')
    plot(xm2y,[zm2y(1),zm2y(1)],'r--')
    plot(xm2y,[zm2y(2),zm2y(2)],'r--')
    plot([xm2y(1),xm2y(1)],zm2y,'r--')
    plot([xm2y(2),xm2y(2)],zm2y,'r--')
hold off
xlabel('x(m)');ylabel('depth(m)');
xlim([0,8]);set(gca,'XTick',(0:8))
daspect([1,1,1]);
ca = caxis;
set(gca,'fontsize',20,'fontname','Times')
set(gca,'Position',[0.05,0.02,0.4,0.5])

i = 2;
f_slice = dir(fullfile(resultdir,sprintf('result_xcorr_????_%02d.???',i-1)));
slice_tag = [sprintf("x=%gm",x(slices_ix(i))),"y(m)","depth(m)"];
[slicex2,slicex_sum2] =  slice_summation(f_slice,y,z,slice_tag,fig_result_dir,slice_reload,fig_para,subplot(2,2,4));
hold on
plot(ym2x,[zm2x(1),zm2x(1)],'r--')
plot(ym2x,[zm2x(2),zm2x(2)],'r--')
plot([ym2x(1),ym2x(1)],zm2x,'r--')
plot([ym2x(2),ym2x(2)],zm2x,'r--')
hold off
xlim([0,5]);set(gca,'XTick',(0:5));
xlabel('y(m)');ylabel('depth(m)');
daspect([1,1,1]);
set(gca,'fontsize',20,'fontname','Times')
set(gca,'Position',[0.5,0.02,0.25,0.5])

slicex = {[]};
slicex_sum = {[]};
i = 1;
f_slice = dir(fullfile(resultdir,sprintf('result_xcorr_????_%02d.???',i-1)));
slice_tag = [sprintf("x=%gm",x(slices_ix(i))),"y(m)","depth(m)"];
[slicex1,slicex_sum1] =  slice_summation(f_slice,y,z,slice_tag,fig_result_dir,slice_reload,fig_para,subplot(2,2,2));
hold on
plot(ym1x,zm1x,'r--')
hold off
ylabel('depth(m)');xlabel('')
xlim([0,5]);set(gca,'XTickLabel','')
daspect([1,1,1]);
set(gca,'fontsize',20,'fontname','Times')
set(gca,'Position',[0.5,0.28,0.25,0.5])

slicez = {[]};
slicez_sum = {[]};
fig_paraz = [-inf, fig_para(2)];
for i =1:length(slices_iz)
    f_slice = dir(fullfile(resultdir,sprintf('result_zcorr_????_%02d.???',i-1)));
    slice_tag = [sprintf("z=%gm",z(slices_iz(i))),"x(m)","y(m)"];
    [slicez{i},slicez_sum{i}] =  slice_summation(f_slice,x,y,slice_tag,fig_result_dir,slice_reload,fig_paraz,subplot(2,2,1));
end
% model
hold on
plot(xm1z,ym1z,'r--',xm2z,ym2z,'r--')
hold off
ylabel('y(m)');xlabel('')
xlim([0,8]);ylim([0,5])
set(gca,'XTickLabel','')
daspect([1,1,1])
set(gca,'fontsize',20,'fontname','Times')
set(gca,'Position',[0.05,0.35,0.4,0.5])
%%
export_fig(fullfile(fig_result_dir,"slices_with_model.png"),'-transparent')
% export_fig(fullfile(fig_result_dir,"slices_with_model.pdf"),'-transparent')

%% Draw wavefields
reload_wvf = true;

amp_para.fac1 = 0;
amp_para.fac2 = 0;
try
    amp_para.srcx = m.srcx;
    amp_para.srcy = m.srcy;
    amp_para.srcz = m.srcz;
catch
    amp_para.fac1 = 0;
    amp_para.fac2 = 0;
    amp_para.srcx = [0];
    amp_para.srcy = [0];
    amp_para.srcz = [0];
end
% [wavefield,wavefield_sum] = wavefield_summation(f_wavefield_corr,fig_result_dir,wx,wy,wz,amp_para,wslicex,wslicey,wslicez,reload_wvf);
% %%
% nfwf = length(f_wavefield_corr);
% i = 1;
% while i <= nfwf
%     ext = f_wavefield_corr(i).name(end-2:end);
%     if ~strcmp(ext,'bin') && ~strcmp(ext,'dat')
%         f_wavefield_corr(i) = [];
%         i = i-1;
%         nfwf = nfwf-1;
%     end
%     i = i+1;
% end
% %% for test wvf reading
% for i=round(length(f_wavefield_corr)/2)
%     fid = fopen(fullfile(f_wavefield_corr(i).folder,f_wavefield_corr(i).name));
%     wvf_tmp = fread(fid,[nz,nx*ny],'float');
%     fclose(fid);
%     wvf_tmp = reshape(wvf_tmp,[nz,ny,nx]);
% %     wvf_sl = reshape(wvf_tmp(49,:),[ny,nx]);
%     wvf_tmp = permute(wvf_tmp,[2,3,1]);
% %     wvfi = load(fullfile(resultdir,f_wavefield_corr(i).name));
% % %     zsl = reshape(wvfi(:,wslicez(1)),wny,wnx);
% %     wvfi = reshape(wvfi,wny,wnx,wnz);
%     imagesc(squeeze(wvf_tmp(:,:,:)));colorbar
%     % wvfi = wavefield{i};
% %     imagesc(squeeze(wvfi(:,:,20))')
%     pause(0.1)
% end

%% test
% for i = 1: length(wavefield)
%     wf = squeeze(wavefield{i}(:,wslicey,:));
%     imagesc(wf');colorbar
%     title(sprintf("%d",i))
%     pause(0.1)
% end

%% gathers
% workdirp2 = 'dm_g0.025m_400MHz_0.5x0.5_2_0.1x0.1_pstd_8';
% workdirf2 = 'dm_g0.025m_400MHz_0.5x0.5_2_0.1x0.1_fdtd_8';
% workdirp4 = 'dm4_g0.025m_400MHz_0.5x0.5_2_0.1x0.1_pstd_8';
% workdirf4 = 'dm4_g0.025m_400MHz_0.5x0.5_2_0.1x0.1_fdtd_8';
% % [pgather2,pgather_std2] = getGather(workdirp2,49,49,41);
% % [fgather2,fgather_std2] = getGather(workdirf2,49,49,41);
% % [pgather4,pgather_std4] = getGather(workdirp4,49,49,41);
% % [fgather4,fgather_std4] = getGather(workdirf4,49,49,41);
% figure(10)
% showGather(pgather2,pgather_std2,fgather2,fgather_std2)
% figure(11)
% showGather(pgather4,pgather_std4,fgather4,fgather_std4)
%%
function [gather,gather_std] = getGather(workdir,nrecx,nrecy,i)
%%
    gatherDir = fullfile(workdir,'Output');
    gatherStdDir = fullfile(workdir,'STD','Output');
    f_gather = dir(fullfile(gatherDir,'merge_gather_*.dat'));
    f_gatherstd = dir(fullfile(gatherStdDir,'merge_gather_*.dat'));
    disp(['Loading ' num2str(i) 'th file:' f_gather(i).name])
    gather1 = load(fullfile(f_gather(i).folder,f_gather(i).name));
    nt_g = size(gather1,2);
    gather2 = reshape(gather1,nrecy,nrecx,nt_g);
    gather = permute(gather2,[2,1,3]);
    
    gather1_std = load(fullfile(f_gatherstd(i).folder,f_gatherstd(i).name));
    gather2_std = reshape(gather1_std,nrecy,nrecx,nt_g);
    gather_std = permute(gather2_std,[2,1,3]);
end
%%
function showGather(gather1,gather1_std,gather2,gather2_std)
%%
%     global dt
%     nrecy = size(gather1,1);
%     nt_g = size(gather1,3);
    subplot(2,2,1)
    xsl = squeeze(gather1(25,:,:));
    imagesc(xsl(:,500:end)')
    colorbar
%     hold on
%     h = plot([2,nrecy-1],[1,1],'--');
    subplot(2,2,3)
    xslStd = squeeze(gather1_std(25,:,:));
    imagesc(xslStd(:,1:end)')
    colorbar
%     hold on
%     hstd = plot([2,nrecy-1],[1,1],'--');
    subplot(2,2,2)
    xsl = squeeze(gather2(25,:,:));
    imagesc(xsl(:,500:end)')
    colorbar
%     hold on
%     h = plot([2,nrecy-1],[1,1],'--');
    subplot(2,2,4)
    xslStd = squeeze(gather2_std(25,:,:));
    imagesc(xslStd(:,1:end)')
    colorbar

    
%     for it = 1:100:nt_g
%         subplot(2,2,1)
%         delete(h)
%         h = plot([2,nrecy-1],[it,it],'r--');
%         subplot(2,2,3)
%         delete(hstd)
%         hstd = plot([2,nrecy-1],[it,it],'r--');
%         subplot(2,2,2)
%         zsl = gather(:,:,it);
%         imagesc(zsl)
%         title(sprintf("t=%gns",dt*it/1e-9));
%         subplot(2,2,4)
%         zslstd = gather_std(:,:,it);
%         imagesc(zslstd)
%         pause(0.1)
%     end
end

%%
function [slice,slice_sum] = slice_summation(f_slice,x1,x2,slice_tag,fig_out_dir,reload,fig_para,figcmd)
    if nargin < 6
        reload = false;
    end
    if nargin < 7
        zmin = -inf;
        enhance = 0;
    else
        zmin = fig_para(1);
        enhance = fig_para(2);
% enhance: figure enhance mode - 0 for no enhance, 1 for agc, 2 for balance
%           positive and negetive values.
    end
    if nargin < 8
        figcmd = false;
    end
    nx1 = length(x1);
    nx2 = length(x2);
    outdir = f_slice(1).folder;
    outf = fullfile(outdir,"slice_sum_" + slice_tag(1) + ".mat");
%     disp(outf)
    if ~exist(outf,'file') || reload
        slice = zeros(nx1,nx2,length(f_slice));
        slice_sum = zeros(nx1,nx2);
        
        nfslice = length(f_slice);
        i = 1;
        while i <= nfslice
            ext = f_slice(i).name(end-2:end);
            if ~strcmp(ext,'bin') && ~strcmp(ext,'dat')
                f_slice(i) = [];
                i = i-1;
                nfslice = nfslice-1;
            end
            i = i+1;
        end
        parfor i = 1:length(f_slice)
            ext = f_slice(i).name(end-2:end);
            disp(['Loading ' num2str(i) 'th file:' f_slice(i).name])
            if strcmp(ext,'dat')
                slicei = load(fullfile(outdir,f_slice(i).name));
            else
                fid = fopen(fullfile(outdir,f_slice(i).name));
                slicei = fread(fid,[nx2,nx1],'float');
                slicei = slicei';
                fclose(fid);
            end
            slice(:,:,i) = slicei;
            slice_sum = slice_sum + slicei;
        end
        save(outf,'slice','slice_sum')
    else
        load(outf)
    end
%     colormap gray
    set(gcf,'Unit','centimeters')
    set(gcf,'Position',[0,0,29.7,21])
    set(gca,'fontsize',30,'fontname','Times')
    if enhance == 1
        slice_sum_en = agc(slice_sum')';
    elseif enhance == 2
        slice_sum_en = slice_sum;
        rmax = max(max(slice_sum));
        rmin = min(min(slice_sum));
        slice_sum_en(slice_sum > 0) = slice_sum(slice_sum > 0)/rmax;
        slice_sum_en(slice_sum < 0) = slice_sum(slice_sum < 0)/-rmin;
    else
        slice_sum_en = slice_sum;
%         colorbar
    end
    
    if figcmd == false
        figure(1)
        clf;
        imagesc(x1,x2(x2>zmin),slice_sum(:,x2>zmin)')
        title("slice at " + slice_tag(1))
        daspect([1,1,1]);colorbar
        xlabel(slice_tag(2));ylabel(slice_tag(3))
        set(gca,'fontsize',24,'fontname','Times')
        export_fig(fullfile(fig_out_dir,"slice_" + slice_tag(1) + "_ori.png"),'-transparent')

        imagesc(x1,x2(x2>zmin),slice_sum_en(:,x2>zmin)')
        daspect([1,1,1]);
        xlabel(slice_tag(2));ylabel(slice_tag(3))
        set(gca,'fontsize',24,'fontname','Times')
        export_fig(fullfile(fig_out_dir,"slice_" + slice_tag(1) + ".png"),'-transparent')
    else
        figure(1)
        imagesc(x1,x2(x2>zmin),slice_sum_en(:,x2>zmin)')
        daspect([1,1,1])
        xlabel(slice_tag(2));ylabel(slice_tag(3))
        set(gca,'fontsize',24,'fontname','Times')
    end
end

%%
function [wavefield,wavefield_sum] = wavefield_summation(f_wavefield_corr,fig_result_dir,x,y,z,amp_para,wslicex,wslicey,wslicez,reload)
    nx = length(x);
    ny = length(y);
    nz = length(z);
    outdir = f_wavefield_corr(1).folder;
    %%
    f_wvf = fullfile(outdir,'result_wavefield.mat');
    f_wvf_sum = fullfile(outdir,'result_wavefield_sum.mat');
    if ~exist(f_wvf_sum,'file') || reload
        if ~exist(f_wvf,'file') || reload
            wavefield = {[]};
            
            nfwf = length(f_wavefield_corr);
            i = 1;
            while i <= nfwf
                ext = f_wavefield_corr(i).name(end-2:end);
                if ~strcmp(ext,'bin') && ~strcmp(ext,'dat')
                    f_wavefield_corr(i) = [];
                    i = i-1;
                    nfwf = nfwf-1;
                end
                i = i+1;
            end
            parfor i = 1:length(f_wavefield_corr)
                ext = f_wavefield_corr(i).name(end-2:end);
                disp(['Loading ' num2str(i) 'th file:' f_wavefield_corr(i).name])
                if strcmp(ext,'dat')
                    wvf_tmp = reshape(load(fullfile(f_wavefield_corr(i).folder,f_wavefield_corr(i).name)),ny,nx,nz);
                    wavefield{i} = permute(wvf_tmp,[2,1,3]);
                else
                    fid = fopen(fullfile(f_wavefield_corr(i).folder,f_wavefield_corr(i).name));
                    wvf_tmp = fread(fid,[nz,ny*nx],'float');
                    fclose(fid);
                    wvf_tmp = reshape(wvf_tmp,[nz,ny,nx]);
                    wavefield{i} = permute(wvf_tmp,[3,2,1]);
                end
            end
            save(f_wvf,'wavefield','-v7.3')
        else
            load(f_wvf)
        end
        %%
        disp('Summing wavefields...')
        for i = 1:length(f_wavefield_corr)
            wavefield_i = wavefield{i};
            if ~(amp_para.fac1==0 && amp_para.fac2==0)
                [wavefield_amp,factor1,factor2] = amp_gain_distance(wavefield_i,[amp_para.srcx(i),amp_para.srcy(i),amp_para.srcz],x,y,z,amp_para.fac1,pi,amp_para.fac2);% mul1&2 are 2 and 5
            end
            if i == 1
                wavefield_sum = wavefield_i;
                if ~(amp_para.fac1==0 && amp_para.fac2==0)
                    wavefield_amp_sum = wavefield_amp;
                    factor1_sum = factor1;
                    factor2_sum = factor2;
                end
            else
                wavefield_sum = wavefield_sum + wavefield_i;
                if ~(amp_para.fac1==0 && amp_para.fac2==0)
                    wavefield_amp_sum = wavefield_amp_sum + wavefield_amp;
                    factor1_sum = factor1_sum + factor1;
                    factor2_sum = factor2_sum + factor2;
                end
            end
        end
        if amp_para.fac1==0 && amp_para.fac2==0
            wavefield_amp_sum = wavefield_sum;
        end
%         wavefield_amp_sum1 = wavefield_amp_sum./factor2_sum;
        if ~(amp_para.fac1==0 && amp_para.fac2==0)
            save(f_wvf_sum,'wavefield_sum','wavefield_amp_sum','factor1_sum','factor2_sum','-v7.3')
        else
            save(f_wvf_sum,'wavefield_sum','wavefield_amp_sum','-v7.3')
        end
    else
        load(f_wvf)
        load(f_wvf_sum)
    end

    %%
    figure(10)
%     colormap gray
    set(gcf,'Unit','centimeters')
    set(gcf,'Position',[0,0,29.7,21])
    set(gca,'fontsize',30,'fontname','Times')
    global xm1z ym1z xm2z ym2z
    if wslicez
        for i = wslicez
            zi = z(i);
            zslice = squeeze(wavefield_sum(:,:,i));
            imagesc(x,y,zslice');
%             title(['zslice at z=' num2str(zi) 'm']);colorbar

            % outline model, should add manully
            hold on
            plot(xm1z,ym1z,'r--',xm2z,ym2z,'r--')
            hold off

            daspect([1,1,1])
            xlabel('x(m)');ylabel('y(m)')
            set(gca,'fontsize',20,'fontname','Times')
            
            colorbar
%             [camin,camax] = caxis;
%             cabdr = (camax-camin)/2;
%             caxis([-cabdr,cabdr]);

            export_fig(fullfile(fig_result_dir,['wvf_zslice at z=' num2str(zi) 'm.png']),'-transparent')
            pause(0.1)
        end
    end
    
    %%
    if wslicex
        for i = wslicex
            xi = x(i);
            xslice = squeeze(wavefield_sum(i,:,:));
            xslice1 = agc(xslice);
            imagesc(y,z,xslice');
%             title(['xslice at x=' num2str(xi) 'm']);colorbar;
            daspect([1,1,1])
            xlabel('y(m)');ylabel('depth(m)')
            set(gca,'fontsize',20,'fontname','Times')

    %         % outline model, should add manully
% r=0.5;
% theta = 0:0.1:2*pi;
% x_t = 2.5+r * cos(theta);
% global zd
% z_t = zd+r * sin(theta);
% hold on
% plot(x_t,z_t,'r--')
% hold off

    colorbar
%     [camin,camax] = caxis;
%     cabdr = (camax-camin)/2;
%     caxis([-cabdr,cabdr]);

            export_fig(fullfile(fig_result_dir,['wvf_xslice at x=' num2str(xi) 'm.png']),'-transparent')
            pause(0.1)
        end
    end
    %%
    if wslicey
        for i = wslicey
            yi = y(i);
            zind = z>-5;
            xind = x>-5;%x>0.4&x<5-0.4;
            yslice = squeeze(wavefield_sum(xind,i,zind));
            op_length = 8;
            yslice1 = agc(yslice);%,1:size(yslice,1),8,op_length/10,1);
            imagesc(x(xind),z(zind),yslice');
%             caxis([-1e-5,1e-5]);
%             title(['yslice at y=' num2str(yi) 'm']);colorbar;

    % % outline model, should add manully
    %         hold on
    %         % draw layer
    %         surf_pos = [0.9,1.5,2.4];
    %         zm1 = ones(size(x))*surf_pos(1);
    %         zm2 = ones(size(x))*surf_pos(2);
    %         zm3 = ones(size(x))*surf_pos(3);
    %         dh = 0.2; %m
    %         for ix = 1:length(x)
    %             xi = x(ix);
    %             z1i = zm1(ix);
    %             z1hi = zm1(ix)+dh;
    %             z2i = zm2(ix);
    %             z2hi = zm2(ix)+dh;
    %             z3i = zm3(ix);
    %             z3hi = zm3(ix)+dh;
    %             dot1i = [xi,yi,z1i];
    %             dot1hi = [xi,yi,z1hi];
    %             dot2i = [xi,yi,z2i];
    %             dot2hi = [xi,yi,z2hi];
    %             dot3i = [xi,yi,z3i];
    %             dot3hi = [xi,yi,z3hi];
    %             
    %             temp1 = dot([7.5,-0.75,-2.5],(dot1i-[2 0 0.8]));
    %             temp1h =  dot([7.5,-0.75,-2.5],(dot1hi-[2 0 0.8]));
    %             if temp1 > 0 && temp1h < 0
    %                 zm1(ix) = (7.5*(xi-2)-0.75*yi)/2.5 + 0.8;
    %             elseif temp1h > 0
    %                 zm1(ix) = z1hi;
    %             end
    %             
    %             temp2 = dot([7.5,-0.75,-2.5],(dot2i-[2 0 0.8]));
    %             temp2h =  dot([7.5,-0.75,-2.5],(dot2hi-[2 0 0.8]));
    %             if temp2 > 0 && temp2h < 0
    %                 zm2(ix) = (7.5*(xi-2)-0.75*yi)/2.5 + 0.8;
    %             elseif temp2h > 0
    %                 zm2(ix) = z2hi;
    %             end
    %             
    %             temp3 = dot([7.5,-0.75,-2.5],(dot3i-[2 0 0.8]));
    %             temp3h =  dot([7.5,-0.75,-2.5],(dot3hi-[2 0 0.8]));
    %             if temp3 > 0 && temp3h < 0
    %                 zm3(ix) = (7.5*(xi-2)-0.75*yi)/2.5 + 0.8;
    %             elseif temp3h > 0
    %                 zm3(ix) = z3hi;
    %             end
    %         end
    %         plot(x,zm1,'r--',x,zm2,'r--',x,zm3,'r--')
    % %         % draw fault: 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0
    % %         zm = (7.5*(x-2)-0.75*yi)/2.5 + 0.8;
    % %         plot(x,zm,'r-.')
    %         plot(x,0.92*ones(size(x)),'k')% show where zslice locate.
    %         hold off
            daspect([1,1,1])
            xlabel('x(m)');ylabel('depth(m)')
            set(gca,'fontsize',20,'fontname','Times')
            
    colorbar
%     [camin,camax] = caxis;
%     cabdr = (camax-camin)/2;
%     caxis([-cabdr,cabdr]);
    
            export_fig(fullfile(fig_result_dir,['wvf_yslice at y=' num2str(yi) 'm.png']),'-transparent')
            pause(0.1)
        end
    end
end
