%%
workdir = 'test1_300MHz_0.6m_0.2m';
resultdir = fullfile(workdir,'Result');
outdir = fullfile(workdir,'.');
f_wavefield_corr = dir(fullfile(resultdir,'result_wavefield_corr.dat*'));

modelfiles = dir(fullfile(workdir,'test*.mat'));
modelfn = fullfile(modelfiles.folder,modelfiles.name);
m = load(modelfn);

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
%%
result_exist = false;
if ~result_exist
    disp(['Loading 1st file...'])
    wavefield_corr_raw = load(fullfile(resultdir,f_wavefield_corr(1).name));
    for i = 2:length(f_wavefield_corr)
        disp(['Loading ' num2str(i) 'th file...'])
        wavefield_corr_raw = wavefield_corr_raw + load(fullfile(resultdir,f_wavefield_corr(i).name));
    end
    wavefield_corr = reshape(wavefield_corr_raw,nz,ny,nx);%z,y,x
    save(fullfile(resultdir,'result_wavefield_corr'),'wavefield_corr','-v7.3')
else
    load(fullfile(resultdir,'result_wavefield_corr'))
end

%%
x = (1:nx)*dx;
y = (1:ny)*dy;
z = ((1:nz)-nz_air)*dz;
%%
for i = 1:nz
    zi = i*dz;
    zslice = squeeze(wavefield_corr(i,:,:));
    imagesc(x,y,zslice);colorbar
    title(['zslice at z=' num2str((i-nz_air)*dz) 'm'])

    % outline model, should add manully
    hold on
    X = x(end);
    Y = y(end);
    Z = z(end);
    r0 = 0.5;
    xm = X/2;
    ym = Y/2;
    zm = 1.5;
    if abs(zi-zm)<=r0
        r = sqrt(r0^2-(zi-zm)^2);
        alpha = 0:0.1:2*pi;
        xx = xm + r * sin(alpha);
        yy = ym + r * cos(alpha);
        plot(xx,yy,'r--')
    end
    hold off
    saveas(gcf,fullfile(outdir,['zslice at z=' num2str((i-nz_air)*dz) 'm.png']))
    pause(0.1)
end
%%
for i = 1:nx
    xi = i*dx;
    xslice = squeeze(wavefield_corr(:,:,i));
    imagesc(y,z,xslice);colorbar
    title(['xslice at x=' num2str(i*dx) 'm'])
    pause(0.1)
    
    % outline model, should add manully
    hold on
    X = x(end);
    Y = y(end);
    Z = z(end);
    r0 = 0.5;
    xm = X/2;
    ym = Y/2;
    zm = 1.5;
    if abs(xi-xm)<=r0
        r = sqrt(r0^2-(xi-xm)^2);
        alpha = 0:0.1:2*pi;
        zz = zm + r * sin(alpha);
        yy = ym + r * cos(alpha);
        plot(yy,zz,'r--')
    end
    hold off
    saveas(gcf,fullfile(outdir,['xslice at x=' num2str(i*dx) 'm.png']))
    pause(0.1)
end

% %% figure zslice
% result_dir = 'Result';
% load(fullfile(result_dir,'result_wavefield_corr'))
% x = dx: dx : (dx * nx);
% y = (1:ny) * dy;
% z = ((1:nz) - nz_air) * dz;
% 
% % z slice
% zslice1 = squeeze(wavefield_corr(nz_air+round(1/nz),:,:));%1m
% zslice2 = squeeze(wavefield_corr(nz_air+round(1.6/nz),:,:));%1.6m
% zslice3 = squeeze(wavefield_corr(nz_air+round(2/nz),:,:));%2m
% figure(1)
% imagesc(y,x,zslice1),title('Depth=1m'),xlabel('y(m)'),ylabel('x(m)'),
% % hold on
% % plot(x,1.2+0 * x,'r--',x,2.4+0 * x,'r--',x,3.6+0 * x,'r--',x,4.8+0 * x,'r--')
% % plot(1.2+0 * x,x,'r--',2.4+0 * x,x,'r--',3.6+0 * x,x,'r--',4.8+0 * x,x,'r--')
% saveas(gca,'zslice_1m.png')
% figure(2)
% imagesc(y,x,zslice2),title('Depth=1.6m'),xlabel('y(m)'),ylabel('x(m)'),hold on
% saveas(gca,'zslice_1.6m.png')
% figure(3)
% imagesc(y,x,zslice3),title('Depth=2m'),xlabel('y(m)'),ylabel('x(m)'),
% % hold on
% % plot(x,1.2+0 * x,'r--',x,2.4+0 * x,'r--',x,3.6+0 * x,'r--',x,4.8+0 * x,'r--')
% % plot(1.2+0 * x,x,'r--',2.4+0 * x,x,'r--',3.6+0 * x,x,'r--',4.8+0 * x,x,'r--')
% saveas(gca,'zslice_2m.png')
% %%
% % x Slice
% xslice1 = squeeze(wavefield_corr(:,:,round(2.4/dx)));%2.4m
% xslice2 = squeeze(wavefield_corr(:,:,round(3/dx)));%3m
% figure(4)
% imagesc(y, z(z > 1),xslice1(z > 1,:)),title('x=2.4m'),xlabel('y(m)'),ylabel('depth(m)'),
% saveas(gca,'xslice_2.4m.png')
% figure(5)
% imagesc(y, z(z > 1),xslice2(z > 1,:)),title('x=3m'),xlabel('y(m)'),ylabel('depth(m)'),
% saveas(gca,'xslice_3m.png')