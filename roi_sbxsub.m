
function plist = roi_sbxsub(img,vimg,lambda,nclust,mag)

% segment a subimage

mimg = mean(img,3);
% th = prctile(mimg(:),50); % minimum average level percentile...

N = size(img,1);
M = size(img,2);

[xx,yy] = meshgrid(1:N,1:M);

img = reshape(img,[N*M size(img,3)]);

img(:,end+1) = lambda*xx(:);
img(:,end+1) = lambda*yy(:);


o = statset;
o.MaxIter = 150;
o.UseParallel = true;

% idx = kmeans(img,nclust,'options',o);
% idx = zeros(size(img,1),6);
% for(nrpt=1:6)
%     idx(:,nrpt) = fkmeans(img,nclust);
% end
% idx = fkmeans(idx,floor(nclust/2));

[uu,~,~] = svd(img,0);
idx = kmeans(uu(:,1:10),nclust);

clf

mimg = (mimg-min(mimg(:)))/(max(mimg(:))-min(mimg(:)))*256;
subplot(1,4,1);
imshow(mimg,gray(256));

subplot(1,4,2)
imagesc(vimg), axis image, axis off;

subplot(1,4,3)
imshow(label2rgb(reshape(idx,[M N])));


%idx = fkmeans(img,nclust);

L = zeros(N,M);
clear plist;

plist = {};

k = 1;
for(i=1:nclust)
%     if length(find(idx == i)) > 20 && length(find(idx == i)) < 500
        bw = reshape(idx==i,[M N]); 
%         subplot(1,4,4),imshow(bw);
        %%
%         Treat original bw differently, to be used for neuropil noisy
%         signal later on.
        %%
        bw = imopen(bw, strel('disk',2,4));
        bw = imclearborder(bw,8);
        bw = imfill(bw,8,'holes');
%         subplot(1,4,4),imshow(bw);
    cc = regionprops(bw,'Area','Solidity','PixelIdxList','Eccentricity','PixelList');
    for(j=1:length(cc))
       mm = mean(mimg(cc(j).PixelIdxList));
       if(cc(j).Area<700*mag^2 && cc(j).Area>50*mag^2 && cc(j).Solidity>0.82 && cc(j).Eccentricity<0.9) % cell size bounds...
            L(cc(j).PixelIdxList)=k;
            plist{k} = cc(j).PixelList;
            k = k+1;
       end
    end
%     end
end

subplot(1,4,4)
imshow(label2rgb(L,'jet','k','shuffle'))

w=1;



