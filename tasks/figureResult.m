%
result_dir = 'Result';
out_dir = '.';
f_wavefield_corr = dir(fullfile(result_dir,'result_wavefield_corr.dat*'));

nx_ori = 100;
ny_ori = 100;
nz_ori = 70;
nz_air_ori = 10;
dx_ori = 0.1;
dy_ori = 0.1;
dz_ori = 0.1;
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
for i = 1:nz
    zslice = squeeze(wavefield_corr(i,:,:));
    imagesc(zslice)
    title(['zslice at iz=' num2str(i)])
    pause(0.1)
end

%% figure zslice
result_dir = 'Result';
load(fullfile(result_dir,'result_wavefield_corr'))
x = dx: dx : (dx * nx);
y = (1:ny) * dy;
z = ((1:nz) - nz_air) * dz;

% z slice
zslice1 = squeeze(wavefield_corr(nz_air+round(1/nz),:,:));%1m
zslice2 = squeeze(wavefield_corr(nz_air+round(1.6/nz),:,:));%1.6m
zslice3 = squeeze(wavefield_corr(nz_air+round(2/nz),:,:));%2m
figure(1)
imagesc(y,x,zslice1),title('Depth=1m'),xlabel('y(m)'),ylabel('x(m)'),
% hold on
% plot(x,1.2+0 * x,'r--',x,2.4+0 * x,'r--',x,3.6+0 * x,'r--',x,4.8+0 * x,'r--')
% plot(1.2+0 * x,x,'r--',2.4+0 * x,x,'r--',3.6+0 * x,x,'r--',4.8+0 * x,x,'r--')
saveas(gca,'zslice_1m.png')
figure(2)
imagesc(y,x,zslice2),title('Depth=1.6m'),xlabel('y(m)'),ylabel('x(m)'),hold on
saveas(gca,'zslice_1.6m.png')
figure(3)
imagesc(y,x,zslice3),title('Depth=2m'),xlabel('y(m)'),ylabel('x(m)'),
% hold on
% plot(x,1.2+0 * x,'r--',x,2.4+0 * x,'r--',x,3.6+0 * x,'r--',x,4.8+0 * x,'r--')
% plot(1.2+0 * x,x,'r--',2.4+0 * x,x,'r--',3.6+0 * x,x,'r--',4.8+0 * x,x,'r--')
saveas(gca,'zslice_2m.png')
%%
% x Slice
xslice1 = squeeze(wavefield_corr(:,:,round(2.4/dx)));%2.4m
xslice2 = squeeze(wavefield_corr(:,:,round(3/dx)));%3m
figure(4)
imagesc(y, z(z > 1),xslice1(z > 1,:)),title('x=2.4m'),xlabel('y(m)'),ylabel('depth(m)'),
saveas(gca,'xslice_2.4m.png')
figure(5)
imagesc(y, z(z > 1),xslice2(z > 1,:)),title('x=3m'),xlabel('y(m)'),ylabel('depth(m)'),
saveas(gca,'xslice_3m.png')