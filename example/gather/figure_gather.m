%%
nx = 250;
ny = 250;

loc = load('merge_gather_loc_0060.dat');
gather = load('merge_gather_0060.dat');
loc_nmo = load('merge_gather_loc_NMO_0060.dat');
gather_nmo = load('merge_gather_NMO_0060.dat');

%%
img = zeros(nx,ny,size(gather,2));
img_nmo = zeros(nx,ny,size(gather_nmo,2));

for i=1:size(loc,1)
    img(loc(i,1),loc(i,2),:) = gather(i,:);
end

for i=1:size(loc_nmo,1)
    img_nmo(loc_nmo(i,1),loc_nmo(i,2),:) = gather_nmo(i,:);
end

%%
for i=70:180
    subplot(1,2,1)
    imagesc(squeeze(img(:,i,:))');colorbar
    title(['image at y\_grid=' num2str(i)])
    subplot(1,2,2)
    imagesc(squeeze(img_nmo(:,i,:))');colorbar
    title(['image after nmo at y\_grid=' num2str(i)])
    pause()
end