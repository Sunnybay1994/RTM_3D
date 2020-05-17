%%
homedir = '3layers_with_0.2m_fault_300MHz_0.4m_0.2m';
workdir = fullfile(homedir,'RTM0_nmo');
indir = fullfile(workdir,'Input');
infn = fullfile(indir,'src.in_0000');
outdir = workdir;
%%
modelfiles = dir(fullfile(homedir,'model.mat'));
modelfn = fullfile(modelfiles.folder,modelfiles.name);
m = load(modelfn);

dx = m.dx;
dy = m.dy;
dz = m.dz;

nx = m.nx;
ny = m.ny;
nz_air = m.nz_air;
nz = m.nz;

dt = m.dt;%ns

%%

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
gather0 = data.data;
%%
field = zeros(nx,ny,nt);
for i=1:nsrc
    field(srcloc(i,1),srcloc(i,2),:) = gather0(i,:);
%     plot(gather0(i,:))
%     title(num2str(i))
%     pause(0.1)
end
%%
for i=ny/2
    imagesc(flipud(squeeze(field(:,i,:))')),colorbar
    title(['iy=' num2str(i)])
%     pause(0.1)
    saveas(gcf,fullfile(outdir,'input.png'));
end