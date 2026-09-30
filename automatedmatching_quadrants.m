function automatedmatching_quadrants(frame_3D4,frame_3D5,savefile,savecond)
% Generates matching points bewteen volumes
% 
% 
% inputs,
%   frame_3D4: structural data of first volume
%   frame_3D5: structural data of second volume
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
scores = zeros(4,1);
% Coordinates for left
[xp,yp] = meshgrid(round((0.225:0.1:0.325)*sx),round((0.2:0.075:0.875)*sy));
scores(1) = calccorr(xp,yp,A,B); 
% Coordinates for right
[xp,yp] = meshgrid(round((0.675:0.1:0.775)*sx),round((0.2:0.075:0.875)*sy)); 
scores(2) = calccorr(xp,yp,A,B); 
% Coordinates for upper 
[xp,yp] = meshgrid(round((0.2:0.075:0.875)*sy), round((0.225:0.1:0.325)*sx));
scores(3) = calccorr(xp,yp,A,B); 
% Coordinates for lower
[xp,yp] = meshgrid(round((0.2:0.075:0.875)*sy), round((0.675:0.1:0.775)*sx));
scores(4) = calccorr(xp,yp,A,B); 

[~,idx] = max(scores);

switch idx
    case 1
        [xp,yp] = meshgrid(round((0.225:0.1:0.325)*sx),round((0.2:0.075:0.875)*sy));
    case 2
        [xp,yp] = meshgrid(round((0.675:0.1:0.775)*sx),round((0.2:0.075:0.875)*sy));
    case 3
        [xp,yp] = meshgrid(round((0.2:0.075:0.875)*sy), round((0.225:0.1:0.325)*sx));
    case 4
        [xp,yp] = meshgrid(round((0.2:0.075:0.875)*sy), round((0.675:0.1:0.775)*sx));
end

locs = [xp(:),yp(:)];
offsetx = round(0.125*sx);
offsety = round(0.125*sy);
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

figure
showMatchedFeatures(A2,B2,locs,matchedpoints,"montag")
title("Matched Inlier Points")

if nargin > 2
    if savecond
        save(savefile,'locs','matchedpoints')
    end
end
end

%% Helper functions
function score = calccorr(xp,yp,A,B)
[sy,sx] = size(A);
offsetx = round(0.125*sx);
offsety = round(0.125*sy);
locs = [xp(:),yp(:)];
score = 0;
for ii = 1:size(locs,1)
    loc = locs(ii,:);
    temp = A(loc(1)-offsety:loc(1)+offsety,loc(2)-offsetx:loc(2)+offsetx);
    c = normxcorr2(temp,B);
    score = score + max(c(:));
end
end