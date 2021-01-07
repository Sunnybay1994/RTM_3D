%% load files and parameters
workdir = '3vs2_800MHz_0.2m_0.04m_fdtd_6';
resultdir = fullfile(workdir,'Result');
outdir = fullfile(workdir,'Result');
f_wavefield_corr = dir(fullfile(resultdir,'result_wavefield_corr*.dat'));

modelfiles = dir(fullfile(workdir,'model.mat'));
modelfn = fullfile(modelfiles.folder,modelfiles.name);
m = load(modelfn);
load(fullfile(workdir,'model_sr.mat'))

dx_ori = m.dx;
dy_ori = m.dy;
dz_ori = m.dz;

nx_ori = m.nx;
ny_ori =m.ny;
nz_air_ori = m.nz_air;
nz_ori = m.nz;

slice_outstop = m.outstep_slice;
x_outstep = m.outstep_x_wavefield;
t_outstep = m.outstep_t_wavefield;

nx = nx_ori / x_outstep;
ny = ny_ori / x_outstep;
nz = nz_ori / x_outstep;
nz_air = nz_air_ori / x_outstep;
dx = dx_ori * x_outstep;
dy = dy_ori * x_outstep;
dz = dz_ori * x_outstep;
x = (1:nx)*dx;
y = (1:ny)*dy;
z = ((1:nz)-nz_air)*dz;

%% behavier
result_exist = false;
amp_exist = false;
draw_xslice = true;
draw_yslice = true;
draw_zslice = true;

%% read middle result of wavefield
out_folder = fullfile(workdir,'Output');
f_wavefield_out = dir(fullfile(out_folder,'wvf_Ey_0004_*.bin'));
for i = 1:length(f_wavefield_out)
    fname = fullfile(out_folder,f_wavefield_out(i).name);
    disp("loading " + fname + "...");
    fo = fopen(fname);
    temp = fread(fo,[nx*ny*nz],'single');
    fclose(fo);
    wvf{i} = reshape(temp,nx,ny,nz);
    imagesc(squeeze(wvf{i}(:,:,12))');colorbar
%     caxis([-0.1,0.1])
%     imagesc(wvf{i})
    pause(0.1)
end

%%
if ~result_exist
    wavefield = {[]};
    parfor i = 1:7%length(f_wavefield_corr)
        disp(['Loading ' num2str(i) 'th file:' f_wavefield_corr(i).name])
        wavefield{i} = reshape(load(fullfile(resultdir,f_wavefield_corr(i).name)),nx,ny,nz);
    end
    save(fullfile(outdir,'result_wavefield'),'wavefield','-v7.3')
else
    load(fullfile(outdir,'result_wavefield'))
end
%%
for i=1:7
wvfi = wavefield{i};
imagesc(squeeze(wvfi(:,40,:))')
pause(0.1)
end
%%
if ~amp_exist
    for i = 1:length(f_wavefield_corr)
        wavefield_i = wavefield{i};
        [wavefield_amp,factor1,factor2] = amp_gain_distance(wavefield_i,[srcx(i),srcy(i),srcz],x,y,z,0,pi,0);% mul1&2 are 2 and 5
        if i == 1
            wavefield_sum = wavefield_i;
            wavefield_amp_sum = wavefield_amp;
            factor1_sum = factor1;
            factor2_sum = factor2;
        else
            wavefield_sum = wavefield_sum + wavefield_i;
            wavefield_amp_sum = wavefield_amp_sum + wavefield_amp;
            factor1_sum = factor1_sum + factor1;
            factor2_sum = factor2_sum + factor2;
        end
    end
    save(fullfile(outdir,'result_wavefield_sum'),'wavefield_sum','wavefield_amp_sum','factor1_sum','factor2_sum','-v7.3')
else
    load(fullfile(outdir,'result_wavefield_sum'))
end
wavefield_amp_sum1 = wavefield_amp_sum./factor2_sum;
%%
if draw_zslice
    for i = round(nz/2);
        zi = (i-nz_air)*dz;
        zslice = squeeze(wavefield_amp_sum(i,:,:));
        imagesc(x,y,zslice);colorbar
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
        
        saveas(gcf,fullfile(outdir,['zslice at z=' num2str(zi) 'm.png']))
        pause(0.1)
    end
end
%%
if draw_xslice
    for i = round(nx/2)
        xi = i*dx;
        xslice = squeeze(wavefield_amp_sum(:,:,i));
        xslice1 = agc(xslice);
        imagesc(y,z,xslice);colorbar
        title(['xslice at x=' num2str(xi) 'm'])

%         % outline model, should add manully
%         hold on
%         % 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0
%         zm = (7.5*(xi-2)-0.75*y)/2.5 + 0.8;
%         plot(y,zm,'r--')
%         hold off
        
        saveas(gcf,fullfile(outdir,['xslice at x=' num2str(xi) 'm.png']))
        pause(0.1)
    end
end
%%
if draw_yslice
    for i = round(ny/2)
        yi = i*dy;
        zind = z>0;
        xind = x>0.4&x<5-0.4;
        yslice = squeeze(wavefield_amp_sum1(zind,i,xind));
        op_length = 8;
        yslice1 = agc(yslice);%,1:size(yslice,1),8,op_length/10,1);
        imagesc(x(xind),z(zind),yslice);colorbar
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
        
        saveas(gcf,fullfile(outdir,['yslice at y=' num2str(yi) 'm.png']))
        pause(0.1)
    end
end

% %% show amp_gain_distance
%     iw=1;
%     wavefield_i = reshape(load(fullfile(resultdir,f_wavefield_corr(iw).name)),nz,ny,nx);
%     alpha = pi/6;
%     wavefield_amp = amp_gain_distance(wavefield_i,[srcx(iw),srcy(iw),srcz],x,y,z);%,2,alpha,5);
%     i = round(ny/2);
%     yi = i*dy;
%     zind = z>0;
%     xind = x>0.4&x<5-0.4;
%     yslice = squeeze(wavefield_amp(zind,i,xind));
%     op_length = 8;
%     yslice1 = agc(yslice,1:size(yslice,1),8,op_length/10,0);
%     imagesc(x(xind),z(zind),yslice);colorbar
% %         caxis([-1e-5,1e-5]);
%     title(['yslice at y=' num2str(yi) 'm'])
%     hold on
%     zz1 = (2.5-x)*cot(alpha);
%     zz2 = (x-2.5)*cot(alpha);
%     plot(x,zz1,'r--',x,zz2,'r--')
%     hold off