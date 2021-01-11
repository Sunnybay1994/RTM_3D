%% load files and parameters
workdir = '2vs3_800MHz_0.2m_0.04m_fdtd_6';
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

slicex = m.slicex;
slicey = m.slicey;
slicez = m.slicez;

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
wslicex = slicex / x_outstep;
wslicey = slicey / x_outstep;
wslicez = slicez / x_outstep;

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

%% Draw 

%% Draw wavefields
reload_wvf = false;
amp_para.fac1 = 0;
amp_para.fac2 = 0;
try
    amp_para.srcx = m.srcx;
    amp_para.srcy = m.srcy;
    amp_para.srcz = m.srcz;
catch
    amp_gain_factor = false;
end
wavefield_summation(f_wavefield_corr,fig_result_dir,wx,wy,wz,amp_gain_factor,wslicex,wslicey,wslicez,reload_wvf)

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
function [wavefield] = wavefield_summation(f_wavefield_corr,fig_result_dir,x,y,z,amp_para,slicex,slicey,slicez,reload_wvf)
    nx = length(x);
    ny = length(y);
    nz = length(z);
    outdir = f_wavefield_corr(1).folder;
    %%
    if ~exist(fullfile(outdir,'result_wavefield_sum.mat'),'file') || reload_wvf
        if ~exist(fullfile(outdir,'result_wavefield.mat'),'file') || reload_wvf
            wavefield = {[]};
            parfor i = 1:length(f_wavefield_corr)
                disp(['Loading ' num2str(i) 'th file:' f_wavefield_corr(i).name])
                wvf_tmp = reshape(load(fullfile(f_wavefield_corr(i).folder,f_wavefield_corr(i).name)),ny,nx,nz);
                wavefield{i} = permute(wvf_tmp,[2,1,3]);
    %             zsl = wavefield{i}(:,:,slicez(1));
    %             imagesc(zsl)
    %             pause(0.1)
            end
            save(fullfile(outdir,'result_wavefield'),'wavefield','-v7.3')
        else
            load(fullfile(outdir,'result_wavefield'))
        end
        %%
        disp('Summing wavefields...')
        for i = 1:length(f_wavefield_corr)
            wavefield_i = wavefield{i};
            if amp_para
                [wavefield_amp,factor1,factor2] = amp_gain_distance(wavefield_i,[amp_para.srcx(i),amp_para.srcy(i),amp_para.srcz],x,y,z,amp_para.fac1,pi,amp_para.fac2);% mul1&2 are 2 and 5
            end
            if i == 1
                wavefield_sum = wavefield_i;
                if amp_para
                    wavefield_amp_sum = wavefield_amp;
                    factor1_sum = factor1;
                    factor2_sum = factor2;
                end
            else
                wavefield_sum = wavefield_sum + wavefield_i;
                if amp_para
                    wavefield_amp_sum = wavefield_amp_sum + wavefield_amp;
                    factor1_sum = factor1_sum + factor1;
                    factor2_sum = factor2_sum + factor2;
                end
            end
        end
        if ~amp_para
            wavefield_amp_sum = wavefield_sum;
        end
%         wavefield_amp_sum1 = wavefield_amp_sum./factor2_sum;
        if amp_para
            save(fullfile(outdir,'result_wavefield_sum'),'wavefield_sum','wavefield_amp_sum','factor1_sum','factor2_sum','-v7.3')
        else
            save(fullfile(outdir,'result_wavefield_sum'),'wavefield_sum','wavefield_amp_sum','-v7.3')
        end
    else
        load(fullfile(outdir,'result_wavefield_sum'))
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
            imagesc(x,y,zslice');colorbar
            title(['zslice at z=' num2str(zi) 'm'])

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

            saveas(gcf,fullfile(fig_result_dir,['wvf_zslice at z=' num2str(zi) 'm.png']))
            pause(0.1)
        end
    end
    
    %%
    if slicex
        for i = slicex
            xi = x(i);
            xslice = squeeze(wavefield_amp_sum(i,:,:));
            xslice1 = agc(xslice);
            imagesc(y,z,xslice');colorbar
            title(['xslice at x=' num2str(xi) 'm'])
            daspect([1,1,1])
            xlabel('y(m)');ylabel('depth(m)')
            set(gca,'fontsize',20,'fontname','Times')

    %         % outline model, should add manully
    %         hold on
    %         % 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0
    %         zm = (7.5*(xi-2)-0.75*y)/2.5 + 0.8;
    %         plot(y,zm,'r--')
    %         hold off

            saveas(gcf,fullfile(fig_result_dir,['wvf_xslice at x=' num2str(xi) 'm.png']))
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
            imagesc(x(xind),z(zind),yslice');colorbar
    %         caxis([-1e-5,1e-5]);
            title(['yslice at y=' num2str(yi) 'm'])

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

            saveas(gcf,fullfile(fig_result_dir,['wvf_yslice at y=' num2str(yi) 'm.png']))
            pause(0.1)
        end
    end
end