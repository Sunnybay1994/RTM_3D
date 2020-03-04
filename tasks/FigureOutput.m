%%
workdir = 'fault_6_6_2';
result_dir = fullfile(workdir,'/RTM/Output');
fxslice = dir(fullfile(result_dir,'xSlice*.dat*'));
fyslice = dir(fullfile(result_dir,'ySlice*.dat*'));
fzslice = dir(fullfile(result_dir,'zSlice*.dat*'));
fwave = dir(fullfile(result_dir,'Wave*.dat*'));
%%
nx=300;ny=240;nz=120;
dt = 3e-11*1e9;%ns
tstep = 5;%wavefiled and xys slice output time step
wstep = 2;%wavefiled output position step

%%
xslice={[]};yslice={[]};zslice={[]};wavefield={[]};
parfor i = 1:length(fyslice)
    xfid = fopen(fullfile(result_dir,fxslice(i).name));
    xcell = textscan(xfid,'%f');
    xslice0 = xcell{1};
    xslice{i} = reshape(xslice0,[nz,ny]);
    fclose(xfid);
    yfid = fopen(fullfile(result_dir,fyslice(i).name));
    ycell = textscan(yfid,'%f');
    yslice0 = ycell{1};
    yslice{i} = reshape(yslice0,[nz,nx]);
    fclose(yfid);
    zfid = fopen(fullfile(result_dir,fzslice(i).name));
    zcell = textscan(zfid,'%f');
    zslice0 = zcell{1};
%     zslice{i} = reshape(zslice0,[ny,nx,2]);
    zslice{i} = reshape(zslice0,[ny,nx]);
    fclose(zfid);
    wfid = fopen(fullfile(result_dir,fwave(i).name));
    wcell = textscan(wfid,'%f');
    wavefield0 = wcell{1};
    wavefield{i} = reshape(wavefield0,[nz/wstep,ny/wstep,nx/wstep]);
    fclose(wfid);
end
save(fullfile(workdir,'result'),'xslice','yslice','zslice','wavefield','-v7.3')
%%
x = (0:wstep:nx)*0.02;
y = (0:wstep:ny)*0.02;
z = ((0:wstep:nz)-20)*0.02;
for i = length(xslice)
    figure(11)
    imagesc(y,z,xslice{i});colorbar;
    hold on
    xm = (nx/2-1)*0.02;
    ym = (0:ny-1)*0.02;
    zm = 2*(1-xm/6-ym/6);
    plot(ym,zm,'r')
    xlabel('y/m');ylabel('z/m');
    ylim([0.1,2])
    title(['xslice t=' num2str(dt*(i-1)*tstep) 'ns'])
    ca = caxis;
    saveas(gca,fullfile(workdir,'xslice.png'))
    
    figure(12)
    imagesc(x,z,yslice{i});colorbar;
    hold on
    xm = (0:nx-1)*0.02;
    ym = (ny/2-1)*0.02;
    zm = 2*(1-xm/6-ym/6);
    plot(xm,zm,'r')
    xlabel('x/m');ylabel('z/m');
    ylim([0.1,2])
    caxis(ca)
    title(['yslice t=' num2str(dt*(i-1)*tstep) 'ns'])
    saveas(gca,fullfile(workdir,'yslice.png'))
    figure(13)
    
%     zslices = zslice{i}(:,:,1);
    zslices = zslice{i}(:,:);
    imagesc(x,y,zslices(:,:));colorbar;
    hold on
    xm = (0:nx-1)*0.02;
    zm = 1-0.02;
    ym = 6*(1-xm/6-zm/2);
    plot(xm,ym,'r')
    xlabel('x/m');ylabel('y/m');
    title(['zslice1 t=' num2str(dt*(i-1)*tstep) 'ns'])
    saveas(gca,fullfile(workdir,'zslice.png'))
%     figure(14)
%     zslices = zslice{i}(:,:,2);
%     imagesc(x,y,zslices(:,:));colorbar;
%     xlabel('x/m');ylabel('y/m');
%     title(['zslice2 t=' num2str(dt*(i-1)*tstep) 'ns'])
%     saveas(gca,fullfile(workdir,'zslice2.png'))
%     pause(0.3)
end 

x=(0:size(wavefield{end},3)-1)*0.04;
y=(0:size(wavefield{end},2)-1)*0.04;
z=((0:size(wavefield{end},1)-1)-10)*0.04;
%%
for i = -3:3
    figure(25+i)
    linewf = squeeze(wavefield{end}(:,round(ny/4+12.5*i),:));
%     linewf = squeeze(wavefield{end}(:,:,i));
    imagesc(x,z,linewf);colorbar;
    xlabel('x(m)');ylabel('z(m)');
    ylim([0.1,2])
    caxis([-1e-6,1e-6])
%     title(i)
    title(['y=' num2str(2.4+0.5*i) 'm'])
    xm = (0:nx/2-1)*0.04;
    ym = 2.4+0.5*i;
    zm = 2*(1-ym/6-xm/6);
    hold on
    plot(xm,zm,'r')
    hold off
    saveas(gcf,fullfile(workdir,['line' num2str(3+i,'%02d') '_3dmig.png']));
%     save(fullfile(workdir,['line' num2str(3+i,'%02d') '_3dmig.mat']),'linewf','x','z');
    pause()
end

%%
result = permute(wavefield{end},[3,2,1]);
figure(21)
slice(x,y,z,permute(result,[2,1,3]),[],[0.9,1.9,2.9,3.9],[]);
shading interp;alpha(0.6);
caxis([-1e-6,1e-6]);colorbar
zlim([-0.4,2]);xlim([0,6]);ylim([0,4.8]);

hold on
[xm,ym] = meshgrid(x,y);
zm = 2*(1-xm/6-ym/6);
zm(zm<0) = 0;
surf(xm,ym,zm);shading interp;alpha(0.3)

xlabel('x(m)');ylabel('y(m)');zlabel('depth(m)');
set(gca,'zdir','reverse','FontSize',24)
hold off
saveas(gcf,fullfile(workdir,'3D_result.fig'));
