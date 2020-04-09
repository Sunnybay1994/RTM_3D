%% load files and parameters
workdir = '.';
result_dir = fullfile(workdir,'/RTM0/Output');
fxslice = dir(fullfile(result_dir,'xSlice*.dat*'));
fyslice = dir(fullfile(result_dir,'ySlice*.dat*'));
fzslice = dir(fullfile(result_dir,'zSlice*.dat*'));
fwave = dir(fullfile(result_dir,'Wave*.dat*'));
outdir = fullfile(workdir,'/RTM0');

modelfiles = dir(fullfile(workdir,'model.mat'));
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
result_exist = true;
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
        saveas(gca,fullfile(outdir,'xslice.png'))

        figure(12)
        imagesc(x,z,agc(yslice{i}));colorbar;
        xlabel('x/m');ylabel('z/m');
        title(['yslice y=2.5m'])
        % outline model, should add manully
        hold on
        yi = 2.5;
        % draw layer
        surf_pos = [0.8,1.3,2.3];
        zm1 = ones(size(x))*0.8;
        zm2 = ones(size(x))*1.3;
        zm3 = ones(size(x))*2.3;
        dh = 0.3; %m
        for i = 1:length(x)
            xi = x(i);
            z1i = zm1(i);
            z2i = zm2(i);
            z3i = zm3(i);
            dot1i = [xi,yi,z1i];
            dot2i = [xi,yi,z2i];
            dot3i = [xi,yi,z3i];
            if dot([7.5,-0.75,-2.5],(dot1i-[2 0 0.8])) > 0
                zm1(i) = zm1(i) + dh; 
            end
            if dot([7.5,-0.75,-2.5],(dot2i-[2 0 0.8])) > 0
                zm2(i) = zm2(i) + dh; 
            end
            if dot([7.5,-0.75,-2.5],(dot3i-[2 0 0.8])) > 0
                zm3(i) = zm3(i) + dh; 
            end
        end
        plot(x,zm1,'r--',x,zm2,'r--',x,zm3,'r--')
        % draw fault: 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0
        zm = (7.5*(x-2)-0.75*yi)/2.5 + 0.8;
        plot(x,zm,'r-.')
        hold off
        saveas(gca,fullfile(outdir,'yslice.png'))

        figure(13)
        imagesc(x,y,agc(zslice{i}')');colorbar;
        xlabel('x/m');ylabel('y/m');
        title(['zslice z=1.2m'])
        % outline model, should add manully
        hold on
        zi = 1.2;
        % 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0
        ym = (7.5*(x-2)-2.5*(zi-0.8))/0.75;
        plot(x,ym,'r--')
        hold off
        saveas(gca,fullfile(outdir,'zslice.png'))
    end
end

%%
x = (1:nx)*dx;
y = (1:ny)*dy;
z = ((1:nz)-nz_air)*dz;
if draw_wavefield_yslice
    figure(15)
    for i = 12:10:113
        yi = i*dy;
        wyslice = squeeze(wavefield{end}(:,i,:));
        wyslice1 = agc(wyslice,1:size(wyslice,1),8);
        imagesc(x,z,wyslice1);colorbar
        xlabel('x/m');ylabel('z/m');
        title(['yslice at y=' num2str(yi) 'm'])

        % outline model, should add manully
        hold on
        % draw layer
        surf_pos = [0.8,1.3,2.3];
        zm1 = ones(size(x))*0.8;
        zm2 = ones(size(x))*1.3;
        zm3 = ones(size(x))*2.3;
        
        % 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0
        dh = 0.3; %m
        for i = 1:length(x)
            xi = x(i);
            z1i = zm1(i);
            z2i = zm2(i);
            z3i = zm3(i);
            dot1i = [xi,yi,z1i];
            dot2i = [xi,yi,z2i];
            dot3i = [xi,yi,z3i];
            if dot([7.5,-0.75,-2.5],(dot1i-[2 0 0.8])) > 0
                zm1(i) = zm1(i) + dh; 
            end
            if dot([7.5,-0.75,-2.5],(dot2i-[2 0 0.8])) > 0
                zm2(i) = zm2(i) + dh; 
            end
            if dot([7.5,-0.75,-2.5],(dot3i-[2 0 0.8])) > 0
                zm3(i) = zm3(i) + dh; 
            end
        end
        plot(x,zm1,'r--',x,zm2,'r--',x,zm3,'r--')
        
        % draw fault
        % 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0
        zm = (7.5*(x-2)-0.75*yi)/2.5 + 0.8;
        plot(x,zm,'r-.')
        hold off
        
        saveas(gcf,fullfile(outdir,['yslice at y=' num2str(yi) 'm.png']))
        pause(0.1)
    end
%     for i = -3:3
%         figure(25+i)
%         linewf = squeeze(wavefield{end}(:,round(ny/4+12.5*i),:));
%     %     linewf = squeeze(wavefield{end}(:,:,i));
%         imagesc(x,z,linewf);colorbar;
%         xlabel('x(m)');ylabel('z(m)');
%         ylim([0.1,2])
%         caxis([-1e-6,1e-6])
%     %     title(i)
%         title(['y=' num2str(2.4+0.5*i) 'm'])
%         xm = (0:nx/2-1)*0.04;
%         ym = 2.4+0.5*i;
%         zm = 2*(1-ym/6-xm/6);
%         hold on
%         plot(xm,zm,'r')
%         hold off
%         saveas(gcf,fullfile(outdir,['line' num2str(3+i,'%02d') '_3dmig.png']));
%         pause(0.1)
%     end
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