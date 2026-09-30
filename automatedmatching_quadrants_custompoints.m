function automatedmatching_quadrants_custompoints(frame_3D4,frame_3D5,coordx,...
    coordy,offset,savefile,savecond)
% Generates matching points bewteen volumes
% 
% 
% inputs,
%   frame_3D4: structural data of first volume
%   frame_3D5: structural data of second volume
%   coordx: Choose normalized x coordinates for matching points of 1st volume 
%   coordy: Choose normalized y coordinates for matching points of 1st volume 
%   offset: Choose how big to make box for normalized cross correlation
%   computation
%   savefile: File to save with matching points between volumes
%   savecond: set true if want matching points saved


%Generate enfaces
A = squeeze(mean(frame_3D4)); A = A/max(A(:));
B = squeeze(mean(frame_3D5)); B = B/max(B(:));
A = vesselness2D(A,1:0.5:2.5, [1;1], 0.5, false);
B = vesselness2D(B,1:0.5:2.5, [1;1], 0.5, false);
A2 = squeeze(mean(frame_3D4)); A2 = A2/max(A2(:));
B2 = squeeze(mean(frame_3D5)); B2 = B2/max(B2(:));

[sy,sx] = size(A);
%1st dimension is horizontal, 2nd is vertical
[xp,yp] = meshgrid(round(coordy*sy), round(coordx*sx));

locs = [xp(:),yp(:)];
offsetx = round(offset*sx);
offsety = round(offset*sy);
matchedpoints = zeros(size(locs));
for ii = 1:size(locs,1)
    loc = locs(ii,:);
    temp = A(loc(1)-offsety:loc(1)+offsety,loc(2)-offsetx:loc(2)+offsetx);
    c = normxcorr2(temp,B);
    [ypeak,xpeak] = find(c==max(c(:)));
    matchedpoints(ii,:) = [ypeak-size(temp,1)/2, xpeak-size(temp,2)/2];
end


temp = locs;
locs(:,1) = locs(:,2); locs(:,2) = temp(:,1);
temp2 = matchedpoints;
matchedpoints(:,2) = matchedpoints(:,1); matchedpoints(:,1) = temp2(:,2);

%Auto filter out badidx
diff = matchedpoints-locs;
id1 = abs(diff(:,1) - median(diff(:,1))) > 10;
id2 = abs(diff(:,2) - median(diff(:,2))) > 10;

% Gets rid of poorly matched points
badidx = id1 | id2; 
if sum(badidx) == length(badidx)
    badidx = [];
end
matchedpoints(badidx,:) = []; locs(badidx,:) = [];


A2 = imadjust(A2);
B2 = imadjust(B2);
figure
showMatchedFeatures(A2,B2,locs,matchedpoints,"montag")
title("Matched Inlier Points")

if savecond
    save(savefile,'locs','matchedpoints')
end
end

