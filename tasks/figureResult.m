%% load files and parameters
workdir = '3layers_with_0.2m_fault_0o_300MHz_0.1x0.5_0_0.1x0.5_fdtd_4_0o';
resultdir = fullfile(workdir,'Result');
f_wavefield_corr = dir(fullfile(resultdir,'result_wavefield_corr*.dat'));

fig_result_dir = workdir;%fullfile(workdir,'Result');

modelfiles = dir(fullfile(workdir,'model.mat'));
modelfn = fullfile(modelfiles.folder,modelfiles.name);
m = load(modelfn);
% load(fullfile(workdir,'model_sr.mat'))

dx = m.dx;
dy = m.dy;
dz = m.dz;

nx = m.nx;
ny =m.ny;
nz_air = m.nz_air;
nz = m.nz;

slices_ix = m.slicex;
slices_iy = m.slicey;
slices_iz = m.slicez;

x = (1:nx)*dx;
y = (1:ny)*dy;
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
wx = (1:wnx)*wdx;
wy = (1:wny)*wdy;
wz = ((1:wnz)-wnz_air)*wdz;
wslicex = slices_ix / x_outstep;
wslicey = slices_iy / x_outstep;
wslicez = slices_iz / x_outstep;

% %% read middle result of wavefield
% out_folder = fullfile(workdir,'RTM','Output');
% f_wavefield_out = dir(fullfile(out_folder,'wvf_Ey_0031_*.bin'));
% for i = 1:length(f_wavefield_out)
%     fname = fullfile(out_folder,f_wavefield_out(i).name);
%     disp("loading " + fname + "...");
%     fo = fopen(fname);
%     temp = fread(fo,[nx*ny*nz],'single');
%     fclose(fo);
%     wvf{i} = reshape(temp,nx,ny,nz);
%     imagesc(squeeze(wvf{i}(:,:,12))');colorbar
% %     caxis([-0.1,0.1])
% %     imagesc(wvf{i})
%     pause(0.1)
% end

%% Draw slices
slice_reload = false;
%%
slicey = {[]};
slicey_sum = {[]};
fig_para = [0.1, true];
for i =1:length(slices_iy)
    f_slice = dir(fullfile(resultdir,sprintf('result_ycorr_????_%02d.dat',i-1)));
    slice_tag = [sprintf("y=%gm",y(slices_iy(i))),"x(m)","depth(m)"];
    [slicey{i},slicey_sum{i}] =  slice_summation(f_slice,x,z,slice_tag,fig_result_dir,slice_reload,fig_para);
end
% model
% dot1 = [1.8 0 1];
% dot2 = [0.2 1.6 0];
% dot3 = [0.8 0 0];
% % dot1 = [1.2 1.6 1];
% % dot2 = [0.2 1.6 0.2];
% % dot3 = [1.8 0 0.2];
% n_plane = cross((dot3-dot1),(dot2-dot1));
% yi = y(slices_iy(i));
% zi = -(n_plane(1)*(x-dot1(1)) + n_plane(2)*(yi-dot1(2)))/n_plane(3)+ dot1(3);
% zi(zi<0.3) = 0.3;
% zi(zi>0.8) = 0.8;
% % zi(zi<0.1) = 0.1;
% % zi(zi>0.9) = 0.9;


% hold on
% plot(x,zi,'r--')
% hold off
% export_fig(fullfile(fig_result_dir,"slice_" + slice_tag(1) + "_with_model.png"))
%%
slicex = {[]};
slicex_sum = {[]};
for i =1:length(slices_ix)
    f_slice = dir(fullfile(resultdir,sprintf('result_xcorr_????_%02d.dat',i-1)));
    slice_tag = [sprintf("x=%gm",x(slices_ix(i))),"y(m)","depth(m)"];
    [slicex{i},slicex_sum{i}] =  slice_summation(f_slice,y,z,slice_tag,fig_result_dir,slice_reload);
end
%%
slicez = {[]};
slicez_sum = {[]};
for i =1:length(slices_iz)
    f_slice = dir(fullfile(resultdir,sprintf('result_zcorr_????_%02d.dat',i-1)));
    slice_tag = [sprintf("z=%gm",z(slices_iz(i))),"x(m)","y(m)"];
    [slicez{i},slicez_sum{i}] =  slice_summation(f_slice,x,y,slice_tag,fig_result_dir,slice_reload);
end

%% Draw wavefields
reload_wvf = false;

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
[wavefield,wavefield_sum] = wavefield_summation(f_wavefield_corr,fig_result_dir,wx,wy,wz,amp_para,wslicex,wslicey,wslicez,reload_wvf);

% %% for test wvf reading
% for i=round(length(f_wavefield_corr)/2)
%     wvfi = load(fullfile(resultdir,f_wavefield_corr(i).name));
% %     zsl = reshape(wvfi(:,wslicez(1)),wny,wnx);
%     wvfi = reshape(wvfi,wny,wnx,wnz);
%     zsl = wvfi(:,:,wslicez(1));
%     imagesc(zsl)
%     % wvfi = wavefield{i};
% %     imagesc(squeeze(wvfi(:,:,20))')
%     pause(0.1)
% end

%%
function [slice,slice_sum] = slice_summation(f_slice,x1,x2,slice_tag,fig_out_dir,reload,fig_para)
    if nargin < 6
        reload = false;
    end
    if nargin < 7
        zmin = -inf;
        ifagc = false;
    else
        zmin = fig_para(1);
        ifagc = fig_para(2);
    end
    nx1 = length(x1);
    nx2 = length(x2);
    outdir = f_slice(1).folder;
    outf = fullfile(outdir,"slice_sum_" + slice_tag(1) + ".mat");
%     disp(outf)
    if ~exist(outf,'file') || reload
        slice = zeros(nx1,nx2,length(f_slice));
        slice_sum = zeros(nx1,nx2);
        parfor i = 1:length(f_slice)
            disp(['Loading ' num2str(i) 'th file:' f_slice(i).name])
            slicei = load(fullfile(outdir,f_slice(i).name))
            slice(:,:,i) = slicei;
            slice_sum = slice_sum + slicei;
        end
        save(outf,'slice','slice_sum')
    else
        load(outf)
    end
    figure(11)
    set(gcf,'Unit','centimeters')
    set(gcf,'Position',[0,0,29.7,21])
    set(gca,'fontsize',30,'fontname','Times')
    if ifagc
        slice_sum = agc(slice_sum')';
    end
    imagesc(x1,x2(x2>zmin),slice_sum(:,x2>zmin)')
%     title("slice at " + slice_tag(1))
    daspect([1,1,1])
    xlabel(slice_tag(2));ylabel(slice_tag(3))
    set(gca,'fontsize',20,'fontname','Times')
    export_fig(fullfile(fig_out_dir,"slice_" + slice_tag(1) + ".png"))
    pause(0.1)
end

%%
function [wavefield,wavefield_sum] = wavefield_summation(f_wavefield_corr,fig_result_dir,x,y,z,amp_para,slicex,slicey,slicez,reload)
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
            parfor i = 1:length(f_wavefield_corr)
                disp(['Loading ' num2str(i) 'th file:' f_wavefield_corr(i).name])
                wvf_tmp = reshape(load(fullfile(f_wavefield_corr(i).folder,f_wavefield_corr(i).name)),ny,nx,nz);
                wavefield{i} = permute(wvf_tmp,[2,1,3]);
    %             zsl = wavefield{i}(:,:,slicez(1));
    %             imagesc(zsl)
    %             pause(0.1)
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
    set(gcf,'Unit','centimeters')
    set(gcf,'Position',[0,0,29.7,21])
    set(gca,'fontsize',30,'fontname','Times')
    if slicez
        for i = slicez
            zi = z(i);
            zslice = squeeze(wavefield_amp_sum(:,:,i));
            imagesc(x,y,zslice');
%             title(['zslice at z=' num2str(zi) 'm']);colorbar

    %         % outline model, should add manully
    %         hold on
    %         % 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0

    %         ym = (7.5*(x-2)-2.5*(zi-0.8))/0.75;
    %         plot(x,ym,'r--')
    %         plot(x,2.52*ones(size(x)),'k')% show where yslice locate.
    %         hold off

            daspect([1,1,1])
            xlabel('x(m)');ylabel('y(m)')
            set(gca,'fontsize',20,'fontname','Times')

            export_fig(fullfile(fig_result_dir,['wvf_zslice at z=' num2str(zi) 'm.png']))
            pause(0.1)
        end
    end
    
    %%
    if slicex
        for i = slicex
            xi = x(i);
            xslice = squeeze(wavefield_amp_sum(i,:,:));
            xslice1 = agc(xslice);
            imagesc(y,z,xslice');
%             title(['xslice at x=' num2str(xi) 'm']);colorbar;
            daspect([1,1,1])
            xlabel('y(m)');ylabel('depth(m)')
            set(gca,'fontsize',20,'fontname','Times')

    %         % outline model, should add manully
    %         hold on
    %         % 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0
    %         zm = (7.5*(xi-2)-0.75*y)/2.5 + 0.8;
    %         plot(y,zm,'r--')
    %         hold off

            export_fig(fullfile(fig_result_dir,['wvf_xslice at x=' num2str(xi) 'm.png']))
            pause(0.1)
        end
    end
    %%
    if slicey
        for i = slicey
            yi = y(i);
            zind = z>0;
            xind = x>0;%x>0.4&x<5-0.4;
            yslice = squeeze(wavefield_amp_sum(xind,i,zind));
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

            export_fig(fullfile(fig_result_dir,['wvf_yslice at y=' num2str(yi) 'm.png']))
            pause(0.1)
        end
    end
end