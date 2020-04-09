%% files and parameters
workdir = 'gly_400MHz';
result_dir = fullfile(workdir,'Output');
outdir = fullfile(workdir,'result');
fxslice = dir(fullfile(result_dir,'xSlice*.dat*'));
fyslice = dir(fullfile(result_dir,'ySlice*.dat*'));
fzslice = dir(fullfile(result_dir,'zSlice*.dat*'));
fwave = dir(fullfile(result_dir,'Wave*.dat*'));

nx0=1220;ny0=220;nz0=130;nz0_air=10;
dx0=0.05;dy0=0.05;dz0=0.05;
dt = 0.0587e-9;
tstep = 4;%wavefiled and xys slice output time step
wstep = 2;%wavefiled output position step
nx=nx0/wstep;ny=ny0/wstep;nz=nz0/wstep;nz_air=nz0_air/wstep;
dx=dx0*wstep;dy=dy0*wstep;dz=dz0*wstep;

%% behavier
result_exist = True;
draw_slices = true;
draw_wavefield_yslice = true;
draw_3D_view = false;

%%
if ~result_exist
    xslice={[]};yslice={[]};zslice={[]};wavefield={[]};
    parfor i = 1:length(fyslice)
        xfid = fopen(fullfile(result_dir,fxslice(i).name));
        xcell = textscan(xfid,'%f');
        xslice0 = xcell{1};
        xslice{i} = reshape(xslice0,[nz0,ny0]);
        fclose(xfid);
        yfid = fopen(fullfile(result_dir,fyslice(i).name));
        ycell = textscan(yfid,'%f');
        yslice0 = ycell{1};
        yslice{i} = reshape(yslice0,[nz0,nx0]);
        fclose(yfid);
        zfid = fopen(fullfile(result_dir,fzslice(i).name));
        zcell = textscan(zfid,'%f');
        zslice0 = zcell{1};
        zslice{i} = reshape(zslice0,[ny0,nx0]);
        fclose(zfid);
        wfid = fopen(fullfile(result_dir,fwave(i).name));
        wcell = textscan(wfid,'%f');
        wavefield0 = wcell{1};
        wavefield{i} = reshape(wavefield0,[nz,ny,nx]);
        fclose(wfid);
    end
    save(fullfile(outdir,'result'),'xslice','yslice','zslice','wavefield','-v7.3')
else
    load(fullfile(outdir,'result'))
end
    
%%
if draw_slices
    x = (1:nx0)*dx0;
    y = (1:ny0)*dy0;
    z = ((1:nz0)-nz0_air)*dz0;
    for i = length(xslice)
        figure(11)
        imagesc(y,z,xslice{i});colorbar;
        xlabel('y/m');ylabel('z/m');
        title(['xslice t=' num2str(dt*(i-1)) 'ns'])
        saveas(gca,fullfile(outdir,'xslice.png'))

        figure(12)
        imagesc(x,z,agc(yslice{i}')');colorbar;
        xlabel('x/m');ylabel('z/m');
        title(['yslice y=5.5m'])
        saveas(gca,fullfile(outdir,'yslice.png'))
        
        figure(13) 
        imagesc(x,y,zslice{i});colorbar;
        xlabel('x/m');ylabel('y/m');
        title(['zslice z=1m'])
        saveas(gca,fullfile(outdir,'zslice.png'))
    end
end

%%
if draw_wavefield_yslice
    x = (1:nx)*dx;
    y = (1:ny)*dy;
    z = ((1:nz)-nz_air)*dz;
    figure(15)
    for i = 1:ny
        yi = i*dy;
        wyslice = squeeze(wavefield{end}(:,i,:));
        wyslice1 = agc(wyslice,1:size(wyslice,1),8);
        imagesc(x,z,wyslice1);colorbar
        title(['yslice at y=' num2str(yi) 'm'])
        saveas(gcf,fullfile(outdir,['yslice at y=' num2str(yi) 'm.png']))
        pause(0.1)
    end
end

%%
info = get(0);
win_size = info.ScreenSize;
zz = z(z>=0);
xx = x(x>=1); 
yy = y(3:17);
result00 = permute(wavefield{end}(:,3:17,:),[1,3,2]);
result00 = result00(z>=0,x>=1,:);
result0 = zeros(length(zz),length(xx),length(y));
parfor lnum =3:17
    h = figure(lnum);
    set(h, 'outerposition', win_size);
    d1 = squeeze(wavefield{end}(:,lnum,:));
    d1 = d1(z>=0,x>=1);
%     plot(d1(:,1))
    imagesc(xx,zz,d1);colorbar
    title(['3D RTM Result at y=' num2str((lnum-1)/2) 'm'],'fontsize',36)
    xlabel('Distance(m)','fontsize',24);ylabel('Depth(m)','fontsize',24)
    saveas(h,[outdir '/line' num2str(lnum,'%02d') '_RTM_result.png']);
% agc (100ns for 200MHz)
    d2 = d1;%agc(d1,zz);
%     imagesc(xx,zz,d4);colorbar
%     title(['Image after AGC at y=' num2str((lnum-1)/2)],'fontsize',36)
%     xlabel('Distance(m)','fontsize',24);ylabel('Depth(m)','fontsize',24)
%     saveas(h,[output_dir '/proc/line' num2str(lnum,'%02d') '_4_Image after AGC.png']);
%Trace equalization
    d2_tr_amp = smooth(mean(abs(d2),1))';
    plot(d2_tr_amp)
    title(['Mean Amplitude each trace at y=' num2str((lnum-1)/2) 'm'],'fontsize',36)
    xlabel('Trace','fontsize',24);ylabel('Amplitude','fontsize',24)
    saveas(h,[outdir '/proc/line' num2str(lnum,'%02d') '_1_AmpTrace.png']);
    d3 = ones(size(d1));
    for i = 1:size(d1,1)
        d3(i,:) = d2(i,:)./d2_tr_amp;
    end
    imagesc(xx,zz,d3);colorbar
    title(['3D RTM Result after Trace Equalization at y=' num2str((lnum-1)/2) 'm'],'fontsize',36)
    xlabel('Distance(m)','fontsize',24);ylabel('Depth(m)','fontsize',24)
    saveas(h,[outdir '/proc/line' num2str(lnum,'%02d') '_1_result after trace equalization.png']);
 
% agc (100ns for 200MHz)
    d4 = agc(d3,zz);
    imagesc(xx,zz,d4);colorbar
    title(['3D RTM Result after AGC at y=' num2str((lnum-1)/2) 'm'],'fontsize',36)
    xlabel('Distance(m)','fontsize',24);ylabel('Depth(m)','fontsize',24)
    saveas(h,[outdir '/proc/line' num2str(lnum,'%02d') '_2_result after AGC.png']);
    savefig(h,[outdir '/proc/line' num2str(lnum,'%02d') '_2_result after AGC']);
    result0(:,:,lnum) = d4;
end
save([outdir '/result0'],'result0')



%% z slice
z0 = 0.9;
imagesc(xx,yy,squeeze(result0(abs(zz-z0)<0.01,:,3:17))')
xlabel('x(m)','fontsize',24);ylabel('y(m)','fontsize',24)
title(['3D RTM Result at z=' num2str(z0) 'm'],'fontsize',36)
set(gca,'ydir','normal','FontSize',20)
caxis([-1,1])
daspect([1 1 1])
zx1=atan((30.5-25)/(6.3-1))*360/2/pi;

%% y slice
y0 = 4;%3,6
yslice = result0(:,:,abs(y-y0)<0.1);
imagesc(xx,zz,yslice)
set(gca,'fontsize',20);
xlabel('x(m)','fontsize',24);ylabel('Depth(m)','fontsize',24)
title(['3D RTM Result at y=' num2str(y0) 'm'],'fontsize',30)
% caxis([-1e7,1e7]);

% for yslice
% daspect([2 1 1]);
% hold on
% plot(x,0.9*ones(size(x)),'--r');
% hold off

% for detail
xlim([22,34]);ylim([0,4]);daspect([1 1 1]);


%% 3D figure
result = permute(result0,[3,2,1]);
figure(21)
slice(xx,y,zz,result,[],[2,4.5,7],[0.9]);
shading interp;alpha(0.6);
% caxis([-1e4,1e4]);colorbar
zlim([0,4]);xlim([10,50]);ylim([1,8]);
set(gca,'fontsize',20);
xlabel('x(m)');ylabel('y(m)');zlabel('depth(m)');
set(gca,'zdir','reverse','FontSize',24)

% %%
% [X,Y,Z] = meshgrid(xx(xx>20&xx<40),y(3:17),zz(zz>0.5&zz<2.5));
% result1 = result(3:17,xx>20&xx<40,zz>0.5&zz<2.5);
% %%
% hold on
% p = patch(isosurface(X,Y,Z,result1,0.2),'FaceAlpha',0.3);
% isonormals(X,Y,Z,result1,p)
% p.FaceColor = 'red';
% p.EdgeColor = 'none';
% % daspect([1 1 1])
% view(3); 
% axis tight
% % camlight 
% % lighting gouraud
% hold off

%%
% hold on
% C = [24,1,1.7;
%     26,1,2.2;
%     24,1.5,1.7;
%     27,1.5,2.2;
% %     24.2,2,1.7;
% %     28.4,2,2.2;
% %     24.4,2.5,1.65;
% %     28.4,2.5,2.2;
% %     24.6,3,1.65
% %     28.4,3,2.2;
% %     24.6,3.5,1.65;
% %     28.4,3.5,2.2;
% %     24.7,4,1.65;
% %     28.6,4,2.2;
% %     25.2,4.5,1.65;
% %     29,4.5,2.2;
% %     25.4,5,1.6;
% %     28.9,5,2.2;
% %     25.4,5.5,1.6;
% %     29.4,5.5,2.2;
% %     26.4,6,1.6;
% %     29.2,6,2.2;
% %     26.4,6.5,1.6;
% %     29.9,6.5,2.3;
% %     26.2,7,1.6;
% %     31.6,7,2.2;
% %     31.6,7.5,2.2;
%     ];
% x1 = C(:,1); y1 = C(:,2); z1 = C(:,3);
% temp = (C'*C)^-1*C'*ones(size(C,1),1);
% a = temp(1);b = temp(2);c = temp(3);
% [xm,ym] = meshgrid(xx(xx>23&xx<35),y);
% zm = (1-a*xm-b*ym)/c;
% xm(:,all(zm>1.6|zm<1.1,1)) = [];
% ym(:,all(zm>1.6|zm<1.1,1)) = [];
% zm(:,all(zm>1.6|zm<1.1,1)) = [];
% logi_y = y1==2|y1==4.5|y1==7;
% plot3(x1(logi_y),y1(logi_y),z1(logi_y),'r.');
% % axis equal
% % xlabel('x(m)');ylabel('y(m)');zlabel('depth(m)');
% % set(gca,'zdir','reverse','FontSize',24)
% % hold on
% surf(xm,ym,zm,'FaceAlpha',0.3);shading interp;
% zx0 = atan(-b/a)*360/2/pi
% qj0 = 180-acos(c/norm([a,b,c]))*360/2/pi
% hold off

%% ricker
f = 100e6;
T = 20e-9;
t0 = 3.7e-9;
t = 0:dt:T;
ricker = (1-2*(pi*f*(t-t0)).^2).*exp(-(pi*f*(t-t0)).^2);
ep = 8;
dzz = dt*3e8/sqrt(ep);
Z = T*3e8/sqrt(ep); 
src = -interp1(0:dzz:Z,ricker,0:dz:Z);
plot(0:dz:Z,src)

%% temp
d1 = squeeze(wavefield{end}(:,3,:));
for i = 1:302
    tr = d1(z<4,i);
%     trr = deconv(src,tr);
%     subplot(2,1,1)
    zzz = zeros(size(tr)); 
    plot(z(z<4),tr)
    hold on
    plot(z(z<4),zzz)
    hold off
%     subplot(2,1,2)
%     plot(trr)
    ylim([-1e7,1e7]);
    pause()
end

