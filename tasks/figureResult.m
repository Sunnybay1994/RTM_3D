%%
result_dir = 'Result';
out_dir = '.';
f_wavefield_corr = dir(fullfile(result_dir,'result_wavefield_corr.dat*'));

nx_ori = 280;
ny_ori = 200;
nz_ori = 130;
nz_air_ori = 10;
dx_ori = 0.05;
dy_ori = 0.05;
dz_ori = 0.05;
x_outstep = 2;
t_outstep = 5;

nx = nx_ori / x_outstep;
ny = ny_ori / x_outstep;
nz = nz_ori / x_outstep;
nz_air = nz_air_ori / x_outstep;
dx = dx_ori * x_outstep;
dy = dy_ori * x_outstep;
dz = dz_ori * x_outstep;
%%
disp(['Loading 1st file...'])
wavefield_corr_raw = load(fullfile(result_dir,f_wavefield_corr(1).name));
for i = 2:length(f_wavefield_corr)
    disp(['Loading ' num2str(i) 'th file...'])
    wavefield_corr_raw = wavefield_corr_raw + load(fullfile(result_dir,f_wavefield_corr(i).name));
end


wavefield_corr = reshape(wavefield_corr_raw,nz,ny,nx);%z,y,x
save(fullfile(result_dir,'result_wavefield_corr'),'wavefield_corr','-v7.3')
%%
result_dir = 'Result';
load(fullfile(result_dir,'result_wavefield_corr'))
x = 0:dx:(nx*dx);
y = 0:dy:(ny*dy);
z = ((0:nz)-nz_air)*dz;
%%
for i = 1:nz
    zslice = squeeze(wavefield_corr(i,:,:));
    imagesc(x,y,zslice)
    title(['zslice at z=' num2str((i-nz_air)*dz) 'm'])
    hold on
    
    zi = (i-nz_air)*dz;
    % outline rectangle
    xm = 3;
    ym = y(end)/2;
    zm = 2.5;
    lx = 3;
    ly = 6;
    lz = 1.2;
    if zi>(zm-lz/2) && zi<(zm+lz/2)
        xx = xm-lx/2:0.1:xm+lx/2;
        plot(xx,(ym-ly/2)*ones(size(xx)),'r--',xx,(ym+ly/2)*ones(size(xx)),'r--')
        yy = ym-ly/2:0.1:ym+ly/2;
        plot((xm-lx/2)*ones(size(yy)),yy,'r--',(xm+lx/2)*ones(size(yy)),yy,'r--')
    end
    % outline sphere 
    r1_0 = 1;
    r2_0 = 0.4;
    xm = 7;
    ym = 3;
    zm = 2.5;
    theta = 0:0.1:2*pi;
    if zi > zm - r1_0 && zi < zm + r1_0
        r1 = sqrt(r1_0^2-(zi-zm)^2);
        xx1 = xm + r1 * sin(theta);
        yy1 = ym + r1 * cos(theta);    
        plot(xx1,yy1,'r--');
        if zi > zm - r2_0 && zi < zm + r2_0
            r2 = sqrt(r2_0^2-(zi-zm)^2);
            xx2 = xm + r2 * sin(theta);
            yy2 = ym + r2 * cos(theta);
            plot(xx2,yy2,'r--');
        end
    end
    % outline cylinder
    r1 = 1;
    r2 = 0.4;
    h = 2;
    xm = 7;
    ym = 7;
    zm = 2.5;
    theta = 0:0.1:2*pi;
    if zi>(zm-h/2) && zi<(zm+h/2)
        xx1 = xm + r1 * sin(theta);
        yy1 = ym + r1 * cos(theta); 
        xx2 = xm + r2 * sin(theta);
        yy2 = ym + r2 * cos(theta);
        plot(xx1,yy1,'r--',xx2,yy2,'r--');
    end
    % outline tetrahedron
    z_top = 1.5;
    h = sqrt(6)*2/3;
    xm = 11;
    ym = 3;
    k_yz = sqrt(3)/(sqrt(6)*2/3);
    k_xy = 2/sqrt(3);
    if zi > z_top && zi < (z_top + h)
        y_max = (zi-z_top)*k_yz;
        yy = (0:0.1:y_max) + ym - 2/3*y_max;
        yy_end = y_max + ym - 2/3*y_max;
        x_max = k_xy * (yy-(ym - 2/3*y_max));
        xx1 = xm - 0.5 * x_max;
        xx2 = xm + 0.5 * x_max;
        xx3 = [xx1 xx2];
        plot(xx1,yy,'r--',xx2,yy,'r--',xx3,yy_end*ones(size(xx3)),'r--')
    end
    ym = 7;
    if zi > z_top && zi < (z_top + h)
        y_max = (z_top + h - zi) * k_yz;
        yy = (0:0.1:y_max) + ym - 2/3*y_max;
        yy_end = y_max + ym - 2/3*y_max;
        x_max = k_xy * (yy-(ym - 2/3*y_max));
        xx1 = xm - 0.5 * x_max;
        xx2 = xm + 0.5 * x_max;
        xx3 = [xx1 xx2];
        plot(xx1,yy,'r--',xx2,yy,'r--',xx3,yy_end*ones(size(xx3)),'r--')
    end
%     saveas(gcf,['zslice at z=' num2str((i-nz_air)*dz) 'm.png'])
    hold off
    pause(0.1)
end
%%
for i = 1:nx
    xslice = squeeze(wavefield_corr(:,:,i));
    imagesc(y,z,xslice)
    title(['xslice at x=' num2str(i*dx) 'm'])
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