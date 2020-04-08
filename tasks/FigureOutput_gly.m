%%
workdir = '400MHz';
result_dir = fullfile(workdir,'Output');
output_dir = fullfile(workdir,'result');
fxslice = dir(fullfile(result_dir,'xSlice*.dat*'));
fyslice = dir(fullfile(result_dir,'ySlice*.dat*'));
fzslice = dir(fullfile(result_dir,'zSlice*.dat*'));
fwave = dir(fullfile(result_dir,'Wave*.dat*'));

%%
nx=1220;ny=219;nz=131;nz_air=12;
dx=0.05;dy=0.05;dz=0.05;
x = (1:nx)*dx;
y = ((1:ny)-1)*dy;
z = ((1:nz)-nz_air)*dz;
dt = 0.0587e-9;
tstep = 4;%wavefiled and xys slice output time step
wstep = 1;%wavefiled output position step
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
    zslice{i} = reshape(zslice0,[ny,nx]);
    fclose(zfid);
    wfid = fopen(fullfile(result_dir,fwave(i).name));
    wcell = textscan(wfid,'%f');
    wavefield0 = wcell{1};
    wavefield{i} = reshape(wavefield0,[nz/wstep,ny/wstep,nx/wstep]);
    fclose(wfid);
end
save(workdir,'xslice','yslice','zslice','wavefield','-v7.3')
%%
for i =1:length(xslice)
    figure(11)
    imagesc(xslice{i});colorbar;
    xlabel('y');ylabel('z');
    title(['t=' num2str(dt*(i-1)*tstep) 'ns'])
    figure(12)
    imagesc(yslice{i});colorbar;
    xlabel('x');ylabel('z');
    title(['t=' num2str(dt*(i-1)*tstep) 'ns'])
    figure(13)
    imagesc(zslice{i});colorbar;
    xlabel('x');ylabel('y');
    title(['t=' num2str(dt*(i-1)*tstep) 'ns'])
    pause(0.1)
end 
%%
% for i =1:length(wavefield)
%     figure(51)
%     linewf = squeeze(wavefield{i}(:,10,:));
%     linewf = linewf(:,200:800);
%     linewf = agc(linewf')';
%     x=1:size(wavefield{end},3)*0.1;
%     z=1:size(wavefield{end},1)*0.1;
%     imagesc(x,z,linewf);colorbar;
%     xlabel('x(m)');ylabel('z(m)');
% %     caxis([-1e-1,1e-1])
%     title(['t=' num2str(i*dt*tstep) 'ns'])
%     pause()
% end
%%
% figure(52)
% for i = 1:184
% sample1=squeeze(wavefield{i}(:,22,:));
% subplot(2,1,1)
% plot(sample1(2,:))
% sample2=pdata{5};
% subplot(2,1,2)
% plot(sample2(end-i*5,:))
% pause(1)
% end

%% just for the ones EW are reversed
% wavefield{end}(:,:,:)=wavefield{end}(:,:,end:-1:1);

%%
% for i =1:40
%     src_cut = 10;
%     x=(1:size(wavefield{end},3))*0.1;
%     x = x(200:800);
%     z=(1:size(wavefield{end},1))*0.1;
%     linewf0 = squeeze(wavefield{end}(src_cut:end,i,:));
%     linewf0 = linewf0(:,200:800);
%     linewf00 = agc(linewf0',z,0.5)';
%     figure(151)    
%     imagesc(x,z(src_cut:end),linewf00(src_cut:end,:));colorbar;
%     xlabel('x(m)');ylabel('z(m)');
% %     caxis([-5e4,5e4])
%     s_y = i*0.2;
%     title(['y=' num2str(s_y)])
%     saveas(gcf,['result/line=' num2str(i,'%02d') '_3dmig.png']);
% %     save(['result/line' num2str(i,'%02d') '_3dmig.mat'],'linewf0','x','z');
% 
%     
%     linewf1_h = smooth(mean(abs(hilbert(linewf0)),2));
%     linewf1 = ones(size(linewf0));
%     for ii = 1:size(linewf0,2)
%         linewf1(:,ii) = linewf0(:,ii)./linewf1_h;
%     end 
%     
%     linewf2_tr_amp = mean(abs(linewf1),1);
% %     plot(linewf1_tr_amp)
%     linewf2 = ones(size(linewf1));
%     for ii = 1:size(linewf1,1)
%         linewf2(ii,:) = linewf1(ii,:);%./linewf2_tr_amp;
%     end
%     
%     linewf = agc(linewf2',z,0.5)';
%     
%     figure(152)
%     imagesc(x,z(src_cut:end),linewf);colorbar;
%     xlabel('x(m)');ylabel('z(m)');
% %     caxis([-1e0,1e0])
%     title(['linenum=' num2str(i)])
%     saveas(gcf,['result/line' num2str(i,'%02d') '_3dmig_p.png']);
% %     save(['result/line' num2str(i,'%02d') '_3dmig_p.mat'],'linewf','x','z');
%     pause(0.1)
% end 

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
    saveas(h,[output_dir '/line' num2str(lnum,'%02d') '_RTM_result.png']);
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
    saveas(h,[output_dir '/proc/line' num2str(lnum,'%02d') '_1_AmpTrace.png']);
    d3 = ones(size(d1));
    for i = 1:size(d1,1)
        d3(i,:) = d2(i,:)./d2_tr_amp;
    end
    imagesc(xx,zz,d3);colorbar
    title(['3D RTM Result after Trace Equalization at y=' num2str((lnum-1)/2) 'm'],'fontsize',36)
    xlabel('Distance(m)','fontsize',24);ylabel('Depth(m)','fontsize',24)
    saveas(h,[output_dir '/proc/line' num2str(lnum,'%02d') '_1_result after trace equalization.png']);
 
% agc (100ns for 200MHz)
    d4 = agc(d3,zz);
    imagesc(xx,zz,d4);colorbar
    title(['3D RTM Result after AGC at y=' num2str((lnum-1)/2) 'm'],'fontsize',36)
    xlabel('Distance(m)','fontsize',24);ylabel('Depth(m)','fontsize',24)
    saveas(h,[output_dir '/proc/line' num2str(lnum,'%02d') '_2_result after AGC.png']);
    savefig(h,[output_dir '/proc/line' num2str(lnum,'%02d') '_2_result after AGC']);
    result0(:,:,lnum) = d4;
end
save([output_dir '/result0'],'result0')



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

