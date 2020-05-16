%% load files and parameters
workdir = '.';
resultdir = fullfile(workdir,'Result');
outdir = fullfile(workdir,'RTM');
f_wavefield_corr = dir(fullfile(resultdir,'result_wavefield_corr.dat*'));

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
draw_xslice = false;
draw_yslice = true;
draw_zslice = true;

%%
if ~result_exist
    for i = 1:length(f_wavefield_corr)
        disp(['Loading ' num2str(i) 'th file:' f_wavefield_corr(i).name])
        wavefield_corr = reshape(load(fullfile(resultdir,f_wavefield_corr(i).name)),nz,ny,nx);
        wavefield_corr_amp = amp_gain_distance(wavefield_corr,[srcx(i),srcy(i),srcz],x,y,z,4);
        if i == 1
            wavefield_corr_sum = wavefield_corr_amp;
        else
            wavefield_corr_sum = wavefield_corr_sum + wavefield_corr_amp;
        end
    end
    save(fullfile(outdir,'result_wavefield_corr'),'wavefield_corr_sum','-v7.3')
else
    load(fullfile(outdir,'result_wavefield_corr'))
end

%%
if draw_zslice
    for i = 1:nz
        zi = (i-nz_air)*dz;
        zslice = squeeze(wavefield_corr_sum(i,:,:));
        imagesc(x,y,zslice);colorbar
        title(['zslice at z=' num2str(zi) 'm'])

        % outline model, should add manully
        hold on
        % 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0
        ym = (7.5*(x-2)-2.5*(zi-0.8))/0.75;
        plot(x,ym,'r--')
        hold off
        
%         saveas(gcf,fullfile(outdir,['zslice at z=' num2str(zi) 'm.png']))
        pause(0.1)
    end
end
%%
if draw_xslice
    for i = 1:nx
        xi = i*dx;
        xslice = squeeze(wavefield_corr_sum(:,:,i));
        xslice1 = agc(xslice);
        imagesc(y,z,xslice1);colorbar
        title(['xslice at x=' num2str(xi) 'm'])
        pause(0.1)

        % outline model, should add manully
        hold on
        % 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0
        zm = (7.5*(xi-2)-0.75*y)/2.5 + 0.8;
        plot(y,zm,'r--')
        hold off
        
%         saveas(gcf,fullfile(outdir,['xslice at x=' num2str(xi) 'm.png']))
        pause(0.1)
    end
end
%%
if draw_yslice
    i=61;
%     wavefield_corr = reshape(load(fullfile(resultdir,f_wavefield_corr(i).name)),nz,ny,nx);
    wavefield_corr_amp = amp_gain_distance(wavefield_corr,[srcx(i),srcy(i),srcz],x,y,z,1,pi/4,2);
    for i = 1:ny%round(ny/2)
        yi = i*dy;
        zind = z>0;
        xind = x>0.4&x<5-0.4;
        yslice = squeeze(wavefield_corr_amp(zind,i,xind));
        op_length = 8;
        yslice1 = agc(yslice,1:size(yslice,1),8,op_length/10,1);
        imagesc(x(xind),z(zind),yslice);colorbar
%         caxis([-1e-5,1e-5]);
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
        
%         % draw fault
%         % 7.5*(x-2)-0.75*y-2.5*(z-0.8)=0
%         zm = (7.5*(x-2)-0.75*yi)/2.5 + 0.8;
%         plot(x,zm,'r-.')
        hold off
        
%         saveas(gcf,fullfile(outdir,['yslice at y=' num2str(yi) 'm.png']))
        pause(0.1)
    end
end