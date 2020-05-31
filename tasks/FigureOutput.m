%% load files and parameters
homedir = 'nmo_vs_full_300MHz_0.4m_0.2m';
workdir = fullfile(homedir,'RTM0');
inputdir = fullfile(workdir,'Input');
result_dir = fullfile(workdir,'Output');
fxslice = dir(fullfile(result_dir,'xSlice*.dat*'));
fyslice = dir(fullfile(result_dir,'ySlice*.dat*'));
fzslice = dir(fullfile(result_dir,'zSlice*.dat*'));
fwave = dir(fullfile(result_dir,'Wave*.dat*'));
outdir = workdir;

modelfiles = dir(fullfile(homedir,'model.mat'));
modelfn = fullfile(modelfiles.folder,modelfiles.name);
m = load(modelfn);

dx_ori = m.dx;
dy_ori = m.dy;
dz_ori = m.dz;

nx_ori = m.nx;
ny_ori =m.ny;
nz_air_ori = m.nz_air;
nz_ori = m.nz;

dt_ori = m.dt;%ns

slicex = m.slicex;
slicey = m.slicey;
slicez = m.slicez;

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

dt = dt_ori/t_outstep;

%% behavier
result_exist = false;
input_exist = false; 
draw_slices = true;
draw_wyslice = false;
draw_wzslice = true;
draw_3D_view = false;

%% figure pre_rtm gathers
if ~input_exist
    fid = fopen(fullfile(inputdir,'src.in_0000'),'r');
    temp = str2num(fgetl(fid));
    ntr = temp(1);nt = temp(2);
    pre_rtm_gathers = zeros(nx_ori,ny_ori,nt);
    
    gather_infos = textscan(fid,'%d %d %d %s\n',ntr);
    srcx = gather_infos{1};
    srcy = gather_infos{2};
    
    for i = 1:ntr
        pre_rtm_gathers(srcx(i),srcy(i),:) = str2num(fgetl(fid));
    end
    fclose(fid);
    save(fullfile(inputdir,'pre_rtm_gathers'),'pre_rtm_gathers','ntr','nt','gather_infos','-v7.3')
else
    load(fullfile(inputdir,'pre_rtm_gathers'))
end
%%
for i = slicey
    yi = i*dy_ori;
    imagesc(flipud(squeeze(pre_rtm_gathers(:,i,:))'))
    title(['y = ' num2str(yi) ' m']);colorbar
%     saveas(gca,fullfile(outdir,['pre_yslice_at_y=' num2str(yi) 'm.png']))
    pause(0.1)
end


%%
if ~result_exist
    xslice={[]};yslice={[]};zslice={[]};wavefield={[]};
    parfor i = 1:length(fyslice)
        xfid = fopen(fullfile(result_dir,fxslice(i).name));
        xcell = textscan(xfid,'%f');
        xslice0 = xcell{1};
        xslice{i} = reshape(xslice0,[nz_ori,ny_ori]);
        fclose(xfid);
        yfid = fopen(fullfile(result_dir,fyslice(i).name));
        ycell = textscan(yfid,'%f');
        yslice0 = ycell{1};
        yslice{i} = reshape(yslice0,[nz_ori,nx_ori]);
        fclose(yfid);
        zfid = fopen(fullfile(result_dir,fzslice(i).name));
        zcell = textscan(zfid,'%f');
        zslice0 = zcell{1};
    %     zslice{i} = reshape(zslice0,[ny,nx,2]);
        zslice{i} = reshape(zslice0,[ny_ori,nx_ori]);
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
    x = (1:nx_ori)*dx_ori;
    y = (1:ny_ori)*dy_ori;
    z = ((1:nz_ori)-nz_air_ori)*dz_ori;
    for i = length(xslice)
        figure(11)
        imagesc(y,z,xslice{i});colorbar;
        xlabel('y/m');ylabel('z/m');
        title(['xslice t=' num2str(dt*(i-1)) 'ns'])
        saveas(gca,fullfile(outdir,['xslice(' num2str(slicex*dx_ori) 'm).png']))

        figure(12)
        imagesc(x,z(z>=0),amp_gain_distance(yslice{i}(z>=0,:),[2.5,2.5,0]));colorbar;
        xlabel('x/m');ylabel('z/m');
%         title(['yslice y=2.5m'])

%         % outline model, should add manully
%         hold on
%         yi = 2.5;
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
%         plot(x,1*ones(size(x)),'k')% show where zslice locate.
%         hold off
        
        daspect([1,1,1])
        xlabel('x(m)');ylabel('depth(m)')
        saveas(gca,fullfile(outdir,['yslice(' num2str(slicey*dy_ori) 'm).png']))

        figure(13)
        imagesc(x,y,zslice{i});colorbar;
        xlabel('x/m');ylabel('y/m');
%         title(['zslice z=1m'])
        
%         % outline model, should add manully
%         hold on
%         zi = 1;
%         % 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0
%         ym = (7.5*(x-2)-2.5*(zi-0.8))/0.75;
%         plot(x,ym,'r--')
%         plot(x,2.5*ones(size(x)),'k')% show where yslice locate.
%         hold off

        daspect([1,1,1])
        xlabel('x(m)');ylabel('y(m)')
        saveas(gca,fullfile(outdir,['zslice(' num2str((slicez-nz_air_ori)*dz_ori) 'm).png']))
    end
end

%%
x = (1:nx)*dx;
y = (1:ny)*dy;
z = ((1:nz)-nz_air)*dz;
%%
if draw_wyslice
    figure(15)
    for i = 12:10:113
        yi = i*dy;
        wyslice = squeeze(wavefield{end}(:,i,:));
        wyslice1 = agc(wyslice);
        imagesc(x,z,wyslice1);colorbar
        xlabel('x/m');ylabel('z/m');
        title(['yslice at y=' num2str(yi) 'm'])

        % outline model, should add manully
        hold on
        % draw layer
        surf_pos = [0.9,1.5,2.4];
        zm1 = ones(size(x))*surf_pos(1);
        zm2 = ones(size(x))*surf_pos(2);
        zm3 = ones(size(x))*surf_pos(3);
        dh = 0.2; %m
        for ix = 1:length(x)
            xi = x(ix);
            z1i = zm1(ix);
            z1hi = zm1(ix)+dh;
            z2i = zm2(ix);
            z2hi = zm2(ix)+dh;
            z3i = zm3(ix);
            z3hi = zm3(ix)+dh;
            dot1i = [xi,yi,z1i];
            dot1hi = [xi,yi,z1hi];
            dot2i = [xi,yi,z2i];
            dot2hi = [xi,yi,z2hi];
            dot3i = [xi,yi,z3i];
            dot3hi = [xi,yi,z3hi];
            
            temp1 = dot([7.5,-0.75,-2.5],(dot1i-[2 0 0.8]));
            temp1h =  dot([7.5,-0.75,-2.5],(dot1hi-[2 0 0.8]));
            if temp1 > 0 && temp1h < 0
                zm1(ix) = (7.5*(xi-2)-0.75*yi)/2.5 + 0.8;
            elseif temp1h > 0
                zm1(ix) = z1hi;
            end
            
            temp2 = dot([7.5,-0.75,-2.5],(dot2i-[2 0 0.8]));
            temp2h =  dot([7.5,-0.75,-2.5],(dot2hi-[2 0 0.8]));
            if temp2 > 0 && temp2h < 0
                zm2(ix) = (7.5*(xi-2)-0.75*yi)/2.5 + 0.8;
            elseif temp2h > 0
                zm2(ix) = z2hi;
            end
            
            temp3 = dot([7.5,-0.75,-2.5],(dot3i-[2 0 0.8]));
            temp3h =  dot([7.5,-0.75,-2.5],(dot3hi-[2 0 0.8]));
            if temp3 > 0 && temp3h < 0
                zm3(ix) = (7.5*(xi-2)-0.75*yi)/2.5 + 0.8;
            elseif temp3h > 0
                zm3(ix) = z3hi;
            end
        end
        plot(x,zm1,'r--',x,zm2,'r--',x,zm3,'r--')
%         % draw fault: 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0
%         zm = (7.5*(x-2)-0.75*yi)/2.5 + 0.8;
%         plot(x,zm,'r-.')
        hold off
        
        saveas(gcf,fullfile(outdir,['yslice at y=' num2str(yi) 'm.png']))
        pause(0.1)
    end
end

%% wv_zslice
if draw_wzslice
    for i = 28
        zi = (i-nz_air)*dz;
        wzslice = squeeze(wavefield{end}(i,:,:));
        imagesc(x,y,wzslice);colorbar
%         title(['zslice at z=' num2str(zi) 'm'])

        % outline model, should add manully
        hold on
        % 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0
        ym = (7.5*(x-2)-2.5*(zi-0.8))/0.75;
        plot(x,ym,'r--')
        plot(x,2.52*ones(size(x)),'k')% show where yslice locate.
        hold off
        daspect([1,1,1])
        xlabel('x(m)');ylabel('y(m)')
        
        saveas(gcf,fullfile(outdir,['zslice at z=' num2str(zi) 'm.png']))
        pause(0.1)
    end
end

%%
if draw_3D_view
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
    saveas(gcf,fullfile(outdir,'3D_result.fig'));
end