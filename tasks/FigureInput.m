%%
infn = fullfile('RTM/Input','src.in_0000');

dx=0.02;
dy=0.02;
dz=0.02;
dt = 3e-11*1e9;%ns

fid = fopen(infn);
info = str2num(fgetl(fid));
nsrc = info(1);
nt = info(2);
srcloc = zeros(nsrc,3);
for i = 1:nsrc
    line = fgetl(fid);
    srcloc(i,:) = str2num(line(1:end-3));
end
fclose(fid);

data = importdata(infn,' ',1+nsrc);
gather = data.data;
%%
for isrc = 1:nsrc
    disp(isrc)
    gather = importdata(['../Output/merge_gather_' num2str(isrc-1,'%04d') '.dat']);
    imagesc(gather');colorbar;
    pause(0.1)
end
%%
x = (0:wstep:nx)*0.02;
y=x;
z = (0:wstep:nz)*0.02;
for i = length(xslice)
    figure(11)
    imagesc(y,z,xslice{i});colorbar;
    xlabel('y/m');ylabel('z/m');
    title(['xslice t=' num2str(dt*(i-1)*tstep) 'ns'])
    saveas(gca,'xslice.png')
    figure(12)
    imagesc(x,z,yslice{i});colorbar;
    xlabel('x/m');ylabel('z/m');
    title(['yslice t=' num2str(dt*(i-1)*tstep) 'ns'])
    saveas(gca,'yslice.png')
    figure(13)
    zslices = zslice{i};
    imagesc(x,y,zslices(:,:));colorbar;
    xlabel('x/m');ylabel('y/m');
    title(['zslice t=' num2str(dt*(i-1)*tstep) 'ns'])
    saveas(gca,'zslice.png')
end 
%%
% for i =1:10
%     figure(15)
%     linewf = squeeze(wavefield{end}(:,-2+5*i,:));
%     x=1:size(wavefield{end},3)*0.1;
%     z=1:size(wavefield{end},1)*0.1;
%     imagesc(x,z,linewf);colorbar;
%     xlabel('x(m)');ylabel('z(m)');
%     caxis([-1e4,1e4])
%     title(['linenum=' num2str(i)])
%     saveas(gcf,['result/line' num2str(i,'%02d') '_3dmig.png']);
%     save(['result/line' num2str(i,'%02d') '_3dmig.mat'],'linewf','x','z');
% %     pause(0.1)
% end 