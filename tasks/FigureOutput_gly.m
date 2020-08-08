%% files and parameters
xspan = 5;yspan = 5;
nxb = 5 + round(xspan/2);nyb = 5 + round(yspan/2); % boundary
workdir_4 = 'gly_400MHz_3.0_0.05_0.1_0.02';
nx0_4 = 1200+10; ny0_4 = 125+10; nz0_4 = 310; nz0_air_4 = 10;
dx0_4 = 0.05; dy0_4 = 0.1; dz0_4 = 0.02;
x0_4 = ((1:nx0_4)-nxb) * dx0_4;
y0_4 = ((1:ny0_4)-nyb) * dy0_4;
z0_4 = ((1:nz0_4) - nz0_air_4) * dz0_4;
dt0_4 = 0.0587e-9;
result_dir_4 = fullfile(workdir_4,'Output');
outdir_4 = fullfile(workdir_4,'result');
mkdir(outdir_4)

workdir_1 = 'gly_100MHz_3.0_0.05_0.1_0.04';
nx0_1 = nx0_4; ny0_1 = 105+10; nz0_1 = 310; nz0_air_1 = 10;
dx0_1 = dx0_4; dy0_1 = dy0_4; dz0_1 = 0.04;
x0_1 = ((1:nx0_1)-nxb) * dx0_1;
y0_1 = ((1:ny0_1)-nyb) * dy0_1;
z0_1 = ((1:nz0_1) - nz0_air_1) * dz0_1;
dt0_1 = 0.1466;
result_dir_1= fullfile(workdir_1,'Output');
outdir_1 = fullfile(workdir_1,'result');
mkdir(outdir_1)


wstep = 2;%wavefiled output position step

nx_4 = round(nx0_4/wstep);
ny_4 = round(ny0_4/wstep);
nz_4 = nz0_4/wstep;
nz_air_4 = nz0_air_4/wstep;
dx_4 = dx0_4 * wstep;
dy_4 = dy0_4 * wstep;
dz_4 = dz0_4 * wstep;
x_4 = ((1:nx_4)-nxb/wstep) * dx_4;
y_4 = ((1:ny_4)-nyb/wstep) * dy_4;
z_4 = ((1:nz_4) - nz_air_4) * dz_4;

nx_1 = round(nx0_1/wstep);
ny_1 = round(ny0_1/wstep);
nz_1 = nz0_1/wstep;
nz_air_1 = nz0_air_1/wstep;
dx_1 = dx0_1 * wstep;
dy_1 = dy0_1 * wstep;
dz_1 = dz0_1 * wstep;
x_1 = ((1:nx_1)-nxb/wstep) * dx_1;
y_1 = ((1:ny_1)-nyb/wstep) * dy_1;
z_1 = ((1:nz_1) - nz_air_1) * dz_1;


% tstep = 4;%wavefiled and xys slice output time step
% dt = dt0 * tstep;



%% load data
result_exist = true;
[wavefields_4,xslice_4,yslice_4,zslice_4] = load_result(result_dir_4,outdir_4,nx0_4,ny0_4,nz0_4,wstep,result_exist);
[wavefields_1,xslice_1,yslice_1,zslice_1] = load_result(result_dir_1,outdir_1,nx0_1,ny0_1,nz0_1,wstep,result_exist);
%% draw slices
% use_agc = true;
% show_all = false;
% draw_slices(outdir_4,xslice_4,yslice_4,zslice_4,x0_4,y0_4,z0_4,use_agc,show_all);
% draw_slices(outdir_1,xslice_1,yslice_1,zslice_1,x0_1,y0_1,z0_1,use_agc,show_all);
%% interp wavefields to look pretty
dxx = 0.05; dyy = dy_1; dzz = 0.04;
[wi4,xx_4,yy_4,zz_4] = interp_wavefield(wavefields_4{end},x_4,y_4,z_4,0,dxx,dyy,dzz,[0.5,x_4(end)-0.5],[],[0,5]);
[wi1,xx_1,yy_1,zz_1] = interp_wavefield(wavefields_1{end},x_1,y_1,z_1,0.8,dxx,dyy,dzz,[0.5,x_1(end)-0.5],[],[0,6]);

%% draw wavefield slices
% draw_wavefield_yslice(outdir_1,wi1,1:round(1/dyy):length(yy_1),xx_1,yy_1,zz_1,false,[])
% draw_wavefield_yslice(outdir_4,wi4,1:round(1/dyy):length(yy_4),xx_4,yy_4,zz_4,false,[])
%%
info = get(0);
win_size = info.ScreenSize;% 获取显示屏的像素尺寸
fh = figure();
set(fh, 'outerposition', win_size);	% 设置图形窗口位置和外尺寸为屏幕大小
%%
z_slice = 0.8;
iz_zslice = find(abs(zz_1-z_slice)<(dzz/2));
draw_wavefield_zslice(outdir_1,wi1,iz_zslice,xx_1,yy_1,zz_1,[27,1,[1,5]],fh);

y_slice = 4;
iy1_yslice = find(abs(yy_1-y_slice)<(dyy/2));
iy4_yslice = find(abs(yy_4-(y_slice+1))<(dyy/2)); % survey lines of 400MHz is 1m south of those of 100MHz
draw_wavefield_yslice(outdir_1,wi1,iy1_yslice,xx_1,yy_1,zz_1,z_slice,[27,1,[6,17]],fh);
draw_wavefield_yslice(outdir_4,wi4,iy4_yslice,xx_4,yy_4,zz_4,false,[27,1,[18,27]],fh);

xlabel('x(m)','fontsize',24);
%%
export_fig(fh,fullfile(outdir_4,['y=' num2str(y_slice) 'm, z=' num2str(z_slice) 'm.png']))


%% 3D figure
if true
    wi1_p = permute(wi1,[2,3,1]);
    figure(21)
    slice(xx_1,yy_1,zz_1,wi1_p,[],[2,5,8],[]);
    daspect([4 1 1])
    axis tight
    shading interp;alpha(0.6);
    % caxis([-1e4,1e4]);colorbar
    zlim([0,4]);xlim([10,50]);ylim([1,9]);
    set(gca,'fontsize',20);
    xlabel('x(m)');ylabel('y(m)');zlabel('depth(m)');
    set(gca,'zdir','reverse','FontSize',24)
    export_fig(fh,fullfile(outdir_1,['3D.png']))
end


%%
function [wavefields,xslice,yslice,zslice] = load_result(result_dir,outdir,nx0,ny0,nz0,wstep,result_exist)
    fprintf('loading result from "%s" to "%s"...\n',result_dir,outdir)
    if ~result_exist
        nx = round(nx0/wstep);
        ny = round(ny0/wstep);
        nz = nz0/wstep;
        fxslice = dir(fullfile(result_dir,'xSlice*.dat*'));
        fyslice = dir(fullfile(result_dir,'ySlice*.dat*'));
        fzslice = dir(fullfile(result_dir,'zSlice*.dat*'));
        fwave = dir(fullfile(result_dir,'Wave*.dat*'));
        xslice={[]};yslice={[]};zslice={[]};wavefields={[]};
        parfor i = 1:length(fyslice)
            fprintf('reading %d/%d...\r',i,length(fyslice))
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
            wavefields0 = wcell{1};
            wavefields{i} = reshape(wavefields0,[nz,ny,nx]);
            fclose(wfid);
        end
        disp('saving...')
        save(fullfile(outdir,'result'),'xslice','yslice','zslice','wavefields','-v7.3')
    else
        disp('loading existing file...')
        r_s = load(fullfile(outdir,'result'));
        xslice = r_s.xslice;
        yslice = r_s.yslice;
        zslice = r_s.zslice;
        wavefields = r_s.wavefields;
    end
end
    
function [fhx,fhy,fhz] = draw_slices(outdir,xslice,yslice,zslice,x0,y0,z0,use_agc,show_all)
    if show_all
        iseq = 1:length(xslice);
    else
        iseq = length(xslice);
    end
    
%     zind = z0>=0;%&z<3;
%     xind = x0>=0;%> 0.5 & x < x(end) - 0.5;
    fhx = figure();
    fhy = figure();
    fhz = figure();
    for i = iseq
        figure(fhx)
        xslicei = xslice{i};
        if use_agc
            xslicei = agc(xslicei);
            fxnote = 'xslice(agc)';
        else
            fxnote = 'xslice';
        end
        imagesc(y0,z0,xslicei);colorbar;
        xlabel('y/m');ylabel('z/m');
        title(fxnote)
        if i == length(xslice)
            saveas(fhx,fullfile(outdir,[fxnote '.png']))
        end

        figure(fhy)
%         imagesc(x0(xind),z0(zind),agc(yslice{i}(zind,xind)));colorbar;
        yslicei = yslice{i};
        if use_agc
            yslicei = agc(yslicei);
            fynote = 'yslice(agc)';
        else
            fynote = 'yslice';
        end
        imagesc(x0,z0,agc(yslicei));colorbar;
        xlabel('x/m');ylabel('z/m');
        title(fynote)
        if i == length(yslice)
            saveas(fhy,fullfile(outdir,[fynote '.png']))
        end
        
        figure(fhz) 
        imagesc(x0,y0,zslice{i});colorbar;
        xlabel('x/m');ylabel('y/m');
        title(['zslice at z=1m'])
        if i == length(xslice)
            saveas(fhz,fullfile(outdir,'zslice.png'))
        end
        
        pause(0.1)
    end
end

function [wavefield_interp,xx,yy,zz] = interp_wavefield(wavefield,x,y,z,z0,dxx,dyy,dzz,xlim,ylim,zlim)
    disp('Interpolating wavefield...')
    [Y,Z,X] = meshgrid(y,z-z0,x);
    xx = x(1):dxx:x(end);
    yy = y(1):dyy:y(end);
    zz = z(1):dzz:z(end);
    if xlim
        xx = xx(xx > xlim(1) & xx < xlim(2));
    end
    if ylim
        yy = yy(yy > ylim(1) & yy < ylim(2));
    end
    if xlim
        zz = zz(zz > zlim(1) & zz < zlim(2));
    end
    [YY,ZZ,XX] = meshgrid(yy,zz,xx);
    wavefield_interp = interp3(Y,Z,X,wavefield,YY,ZZ,XX);
end

function fh = draw_wavefield_zslice(outdir,wavefield,izseq,xx,yy,zz,subp_para,figh)
    if ~isempty(subp_para)
        fh = figh;
        subplot(subp_para(1),subp_para(2),subp_para(3:end))
    else
        info = get(0);
        win_size = info.ScreenSize;% 获取显示屏的像素尺寸
        fh = figure();
        set(fh, 'outerposition', win_size);	% 设置图形窗口位置和外尺寸为屏幕大小
    end
    
    for iz = izseq
        imagesc(xx,yy,squeeze(wavefield(iz,:,:)))
        if isempty(subp_para)
            xlabel('x(m)','fontsize',24);
        end
        ylabel('y(m)','fontsize',24)
    %     title(['3D RTM Result at z=' num2str(z_slice) 'm'],'fontsize',36)
        set(gca,'ydir','normal','FontSize',20)
        daspect([1 1 1])
        if isempty(subp_para)
            export_fig(fh,fullfile(outdir,['zslice at z=' num2str(zz(iz)) 'm.png']))
            pause(0.1)
        else
%             set(gca,'XTick',[])
%             set(gca,'LooseInset',get(gca,'TightInset'))
        end
    end
    % % zx1=atan((30.5-25)/(6.3-1))*360/2/pi;% 100MHz
    % % zx2=atan((31-26)/(9.4-1))*360/2/pi;% 400MHz
end

function fh = draw_wavefield_yslice(outdir,wavefield,iyseq,xx,yy,zz,aux_line,subp_para,figh)
    if ~isempty(subp_para)
        fh = figh;
        subplot(subp_para(1),subp_para(2),subp_para(3:end))
    else
        info = get(0);
        win_size = info.ScreenSize;% 获取显示屏的像素尺寸
        fh = figure();
        set(fh, 'outerposition', win_size);	% 设置图形窗口位置和外尺寸为屏幕大小
    end
    
    for i = iyseq
%         iline = yy(i)-yy(1);
        wyslice = squeeze(wavefield(:,i,:));
%         wyslice1 = zscore(wyslice')';
%         wyslice1 = trace_equal(wyslice1);
        op_length = zz(end)/8;
        [wyslice2,envsm,ma] = agc(wyslice,zz,op_length,op_length/10,0);
        imagesc(xx,zz,wyslice2);
%         colorbar;
        set(gca,'fontsize',20);
        if isempty(subp_para)
            daspect([4 1 1]);
            xlabel('x(m)','fontsize',24);
        else
            daspect([4 1 1]);
%             set(gca,'XTick',[])
%             set(gca,'LooseInset',get(gca,'TightInset'))
        end
        ylabel('Depth(m)','fontsize',24)
%         title(['3D RTM Result at line' num2str(iline) '(y=' num2str(yy(i)) 'm)'],'fontsize',30)
        
        % for yslice
        if aux_line
            hold on
            plot(xx,aux_line*ones(size(xx)),'--r','LineWidth',2);
            hold off
        end
        if isempty(subp_para)
            export_fig(gcf,fullfile(outdir,['yslice at y=' num2str(yy(i)) 'm.png']))
            pause(0.1)
        end
        
%         % for detail
%         xlim([22,34]);
%         daspect([1 1 1]);
%         saveas(gcf,fullfile(outdir,['detailed yslice at y=' num2str(yy(i)) 'm.png']))
        
    end
end

% 
% %%
% [X,Y,Z] = meshgrid(xx(xx>20&xx<40),y(3:17),zz(zz>0.5&zz<2.5));
% result1 = result(3:17,xx>20&xx<40,zz>0.5&zz<2.5);
% %%
% hold on
% p = patch(isosurface(X,Y,Z,result1,0.2),'FaceAlpha',0.3);
% isonormals(X,Y,Z,result1,p)
% p.FaceColor = 'red';
% p.EdgeColor = 'none';
% daspect([1 1 1])
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
% 
% %% ricker
% f = 100e6;
% T = 20e-9;
% t0 = 3.7e-9;
% t = 0:dt:T;
% ricker = (1-2*(pi*f*(t-t0)).^2).*exp(-(pi*f*(t-t0)).^2);
% ep = 8;
% dzz = dt*3e8/sqrt(ep);
% Z = T*3e8/sqrt(ep); 
% src = -interp1(0:dzz:Z,ricker,0:dz:Z);
% plot(0:dz:Z,src)
% 
% %% temp
% d1 = squeeze(wavefield{end}(:,3,:));
% for i = 1:302
%     tr = d1(z<4,i);
% %     trr = deconv(src,tr);
% %     subplot(2,1,1)
%     zzz = zeros(size(tr)); 
%     plot(z(z<4),tr)
%     hold on
%     plot(z(z<4),zzz)
%     hold off
% %     subplot(2,1,2)
% %     plot(trr)
%     ylim([-1e7,1e7]);
%     pause()
% end

