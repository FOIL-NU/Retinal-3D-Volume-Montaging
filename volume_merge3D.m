function [F, F2] = volume_merge3D(homedir,framenum,coordfiles,startcoord,endcoord,...
    deformable,refvol,dp,sx)
% Merges multiple 3D OCT volumes
% 
% 
% inputs,
%   homedir: Folder name with all files
%   framenum: Folder name with each individual single volume
%   coordfiles: File with matching points between volumes
%   startcood: 1st of overlapping volumes
%   endcoord: overlapping volume matched with startcoord
%   deformable: true if desiring deformable registration. false is not
%   refvol: Reference volume whose reference frame to use
%   dp: pixel sizes of each single volume
%   sx: pixel sizes of reconstructed volumes
%
% outputs,
%   F: 3D reconstructed volume using max projections
%   F2: 3D reconstructed volume using mean averaging



%% Load the volumes 
numvols = length(framenum);
frames = cell(1,numvols);

for ii = 1:numvols
    load([homedir num2str(framenum(ii)) '\frame_3D.mat'],'frame_3D')
    frames{ii} = frame_3D;
end
clear frame_3D

%% Load segmentation files
segs = cell(1,numvols);
for ii = 1:numvols
    segs{ii} = loadseg([homedir num2str(framenum(ii))]);
end

%% Create point clouds 
surfs_ilm = cell(1,numvols);
surfs_rpe = cell(1,numvols);
surfs_full = cell(1,numvols);
lmat = cell(1,numvols);

for ii = 1:numvols
    %RPE point cloud
    l1 = medfilt2(depthprofile(segs{ii} == 125)); l1(l1<1) = NaN;
    l1 = refinesurf(l1);
    lmat{ii} = l1;
    surfs_rpe{ii} = findsurf(l1,dp(1,:));
    %ILM point cloud
    l2 = medfilt2(depthprofile(segs{ii} == 255)); l2(l2<1) = NaN;
    l2 = refinesurf(l2);
    surfs_ilm{ii} = findsurf(l2,dp(1,:));
    %Combined point cloud
    surfs_full{ii} = [findsurf(l1,dp(1,:)); findsurf(l2,dp(1,:))];
end
clear segs l1 l2 %Clear unused variables

%% Let's align each surface
t = cell(1,length(coordfiles));
for ii = 1:length(coordfiles)
    t{ii} = pcalignrigid([homedir char(coordfiles(ii))],...
        lmat{startcoord(ii)},lmat{endcoord(ii)},dp(1,:),dp(2,:),...
        surfs_full{startcoord(ii)},surfs_ilm{endcoord(ii)},surfs_rpe{endcoord(ii)});
end

close all

%% Input 
% Map coordinates until they are all in a common reference frame (volume 1)
%%%%%%%%% Will need to change tref to match own data. This are the
%%%%%%%%% transformations required to map to reference volume (volume 1).
%%%%%%%%% For example, 6th volume (644) maps to 5th volume (640), which
%%%%%%%%% maps to 1st volume (632). Hence tref{6} =
%%%%%%%%% rigidtform3d(t{4}.A*t{5}.A); t{5} maps 6th to 5th volume and t{4}
%%%%%%%%% the 5th volume to first volume. 
% Automatic version with this implemented with graph theory coming out.
% Does not need to be changed for demo data.
[maporder,addorder] = mergeorder(startcoord,endcoord,refvol);

tref = cell(1,numvols);
for ii = 1:numvols
    tref{ii} = rigidtform3d(eye(4));
    if ~isempty(maporder{ii})
        order = maporder{ii};
        for jj = length(order):-1:1
            if order(jj) > 0
                tref{ii} = rigidtform3d(t{order(jj)}.A *tref{ii}.A);
            else
                tref{ii} = rigidtform3d(inv(t{abs(order(jj))}.A) *tref{ii}.A);
            end
        end
    end
end

%%%%%%%%

% Show culmative point cloud
surffull_pc = cell(1,numvols);

inlierpoints = [];
for ii = 1:numvols
    surffull_pc{ii} = transformPointsForward(tref{ii}, surfs_full{ii});
    inlierpoints = [inlierpoints;surffull_pc{ii}]; %Update combined pc
end
pcshow(pointCloud(inlierpoints(1:8:end,:)))

%% Determine total boundaries for the volume
dims = size(permute(frames{1},[2 3 1]));
boxes = cell(1,numvols);
compound = [];
for ii = 1:numvols
    [~,box] = pcbox(dims,dp(1,:));
    boxes{ii} = transformPointsForward(tref{ii},box);
    compound = [compound;boxes{ii}];
end

xranges = [min(compound(:,1)), max(compound(:,1))];
yranges = [min(compound(:,2)), max(compound(:,2))];
zranges = [min(compound(:,3)), max(compound(:,3))];

%% Let's try to merge volumes and see how much we are off
xpix = round((xranges(2) - xranges(1))/sx(1));
ypix = round((yranges(2) - yranges(1))/sx(2));
zpix = round((zranges(2) - zranges(1))/sx(3));

% Original reference frame
Rout = imref3d([xpix, ypix, zpix],xranges,yranges,zranges);
%Adding additional frames
frame_3D = permute(frames{1},[2 3 1]);
[F,~] = imwarp(frame_3D,referencegen(size(frame_3D),dp(1,:)),tref{1},'OutputView',Rout);
G = F > 0;
F2 = double(F);
%Translate original image
for ii = 2:numvols
    frame_3D = permute(frames{ii},[2 3 1]);
    [C,~] = imwarp(frame_3D,referencegen(size(frame_3D),dp(2,:)),tref{ii},'OutputView',Rout);
    G1 = C > 0;
    F = blendRegions(F,C);
    G = G + G1;
    F2 = F2 + double(C);
end
F2 = F2 ./ (G+ 0.001);

%% Deformable registration
if deformable
    %Initiate variables
    D = cell(1,numvols);
    ref_ilm = cell(1,numvols);
    ref_rpe = cell(1,numvols);
    rot_ilm = cell(1,numvols);
    rot_rpe = cell(1,numvols);

    %Initiate initial PC
    totcloud_ilm = surfs_ilm{refvol};
    totcloud_rpe = surfs_rpe{refvol};
    ref_ilm{refvol} = surfs_ilm{refvol};
    ref_rpe{refvol} = surfs_rpe{refvol};
    
    for ii = 1:length(addorder)
        ind = addorder(ii);
        [D{ind},ref_ilm{ind},ref_rpe{ind},rot_ilm{ind},rot_rpe{ind}] = calcdef4_2layers(...
            surfs_ilm{ind}, surfs_rpe{ind}, tref{ind}, totcloud_ilm, totcloud_rpe,...
            [400,400,400],xranges,yranges,zranges);

        totcloud_ilm = [totcloud_ilm;ref_ilm{ind}];
        totcloud_rpe = [totcloud_rpe; ref_rpe{ind}];
    end
    pcshow(totcloud_ilm)

    Rout = imref3d([xpix, ypix, zpix],xranges,yranges,zranges);
    Dtemp = zeros(xpix, ypix, zpix,3);
    frame_3D = permute(frames{refvol},[2 3 1]);
    [F,~] = imwarp(frame_3D,referencegen(size(frame_3D),dp(1,:)),tref{refvol},'OutputView',Rout);
    G = F > 0;
    F2 = double(F);
    for ii = 1:length(addorder)
        ind = addorder(ii);
        frame_3D = permute(frames{addorder(ii)},[2 3 1]);
        [C,~] = imwarp(frame_3D,referencegen(size(frame_3D),dp(2,:)),...
            tref{addorder(ii)},'OutputView',Rout);
        
        D2 = D{addorder(ii)};
        for zz = 1:3
            Dtemp(:,:,:,zz) = imresize3(D2(:,:,:,zz),[xpix, ypix, zpix]);
        end
        C = imwarp(C,Dtemp);
        G1 = C > 0;
        F = blendRegions(F,C);
        G = G + G1;
        F2 = F2 + double(C);
    end
    F2 = F2 ./ (G+ 0.001);
end
end


%% Let's run the helper functions
function l = refinesurf(l)
BW = l>0; BW = bwareafilt(BW,1);
l(~BW) = NaN;
end

function [maporder,reforder] = mergeorder(start,final,refframe)
% Construct our graph between nodes
G = graph(start,final);
d = zeros(length(start)+1,1);
% Construct our sequence to reference frame
maporder = cell(1,length(start)+1);

for ii = 1:length(start)+1
    % Shortest path to reference frame
    [P, d(ii)] = shortestpath(G, ii, refframe);
    if length(P) > 1
        temp = zeros(1,d(ii));
        count = 0;
        for jj = length(P)-1:-1:1
            count = count + 1;
            [~,idx] = find(start == P(jj) & final == P(jj+1));
            if ~isempty(idx)
               temp(count) = idx; 
            else
                [~,idx] = find(start == P(jj+1) & final == P(jj));
                temp(count) = -idx; 
            end            
        end
        maporder{ii} = temp;
    end
end

% Progressively add volumes for deformable registration
[~,order] = sort(d);
reforder = order(2:end);
end 


function seg = loadseg(folder)
%Function geared towards loading segmentation data
%Setting home directory
homedir = pwd;
cd(folder)
filenames = dir('*mask*');
%Loading segmentation used
for ii = 1:length(filenames)
    number_str = sprintf('%03d', ii);
    temp = dir(['bscan_' number_str '_mask.png']);
    seg(:,:,ii) = imread(temp.name);
end
%Returning to home directory
cd(homedir)
end

function[firstind,lastind] = depthprofile(BW)
%Top and bottom of our segmentation
firstind = zeros(size(BW,2),size(BW,3));
lastind = firstind;
for ii = 1:size(BW,2)
    for jj = 1:size(BW,3)
        if sum(BW(:,ii,jj)) > 0
            firstind(ii,jj) = find(BW(:,ii,jj),1,'first');
            %lastind(ii,jj) = find(BW(:,ii,jj),1,'last');
        end
    end
end
end


function [D,revcloud_ilm,revcloud_rpe,rot_ilm,rot_rpe] = calcdef4_2layers(surfpoint_ilm,surfpoint_rpe,tform,...
    combinedPC_ilm,combinedPC_rpe,sz,xrange,yrange,zrange)
rot_ilm = zeros(1,size(surfpoint_ilm,1));
rot_rpe = zeros(1,size(surfpoint_rpe,1));

%Let's create a combined PC
surfpoint = [surfpoint_ilm; surfpoint_rpe];
combinedPC = [combinedPC_ilm; combinedPC_rpe];
%Let's run perturbation theory on the combined PC
box = pcboxgen([sz(1) sz(2) sz(3)],[xrange(1) xrange(2)],[yrange(1) yrange(2)],[zrange(1) zrange(2)]);
boxcoord = transformPointsInverse(tform,box);
box2 = box;

%Find the overlapping region in the moving volume
transsurf = transformPointsForward(tform, surfpoint);
%Get rid of overlapping exact coordinates 
[~,~,ib] = intersect(transsurf,combinedPC,"rows");
combinedPC(ib,:) = [];
%Randomly take 10% of points to speed up step
sample = randperm(size(combinedPC,1),round(size(combinedPC,1)/10));
tempcomb = combinedPC(sample,:);
k = boundary(tempcomb(:,1),tempcomb(:,2),1);
in = inpolygon(transsurf(:,1),transsurf(:,2),tempcomb(k,1),tempcomb(k,2));
transsurf = transsurf(in,:);
temp = transformPointsInverse(tform, transsurf);

% %Overlapping region in fixed volume
tempsurf = transsurf; %transsurf(1:2:end,:);
refpoints2 = (tempsurf - mean(tempsurf))*1.3 + mean(tempsurf);
k = boundary(refpoints2(:,1),refpoints2(:,2),1);
in = inpolygon(combinedPC_ilm(:,1),combinedPC_ilm(:,2),refpoints2(k,1),refpoints2(k,2));
combinedPC_ilm = combinedPC_ilm(in,:);
in = inpolygon(combinedPC_rpe(:,1),combinedPC_rpe(:,2),refpoints2(k,1),refpoints2(k,2));
combinedPC_rpe = combinedPC_rpe(in,:);


% Align the subvolumes
y = unique(surfpoint(:,1)); 
y2 = unique(round(temp(:,1)));

% Maxium and minimum indices
minii = find(round(y) == min(y2));
maxii = find(round(y) == max(y2));

if nargout > 1
    revcloud_ilm = surfpoint_ilm;
    revcloud_rpe = surfpoint_rpe;
end

for ii = minii+4:2:maxii-2
    ii
    if ii == minii+4
        lowy = -inf; %y(1) - 260; 
        highy = y(ii) + 260;
      elseif ii == maxii
        lowy = y(ii) - 260;
        highy = inf; %max(y) + 260;
    else
        lowy = y(ii) - 260; 
        highy = y(ii) + 260;
    end
    idx = temp(:,1) > lowy & temp(:,1) < highy;
    if sum(idx) > 0
        if ii == minii+4
            lowy = -inf; %y(1) - 260;
            highy = y(ii) + 130;
        elseif ii == maxii
            lowy = y(ii) - 130;
            highy = inf; %max(y) + 260;
        else
            lowy = y(ii) - 130;
            highy = inf; %y(ii) + 130;
        end
        idx2 = boxcoord(:,1) > lowy & boxcoord(:,1) < highy;
        
        subcloud = temp(idx,:);
        [~,temptform] = genpartialcloud2(subcloud,tform,combinedPC_ilm,...
            combinedPC_rpe);

        ang  = rotm2axang(temptform.R);
        deg = rad2deg(ang(4));

        box2(idx2,:) = transformPointsForward(temptform, boxcoord(idx2,:));
        if nargout > 1
            idx = surfpoint_ilm(:,1) > lowy & surfpoint_ilm(:,1) < highy;
            revcloud_ilm(idx,:) = transformPointsForward(temptform, surfpoint_ilm(idx,:));
            rot_ilm(idx) = deg;

            idx = surfpoint_rpe(:,1) > lowy & surfpoint_rpe(:,1) < highy;
            revcloud_rpe(idx,:) = transformPointsForward(temptform, surfpoint_rpe(idx,:));
            rot_rpe(idx) = deg;
        end
    end
end

%Alterations in each of the coordinates
defx = box(:,1) - box2(:,1);
defy = box(:,2) - box2(:,2);
defz = box(:,3) - box2(:,3);
D = zeros(sz(1), sz(2), sz(3),3); %Create deformation field
D(:,:,:,1) = reshape(defx,sz(1), sz(2), sz(3))/(xrange(2)-xrange(1))*sz(2);
D(:,:,:,2) = reshape(defy,sz(1), sz(2), sz(3))/(yrange(2)-yrange(1))*sz(1);
D(:,:,:,3) = reshape(defz,sz(1), sz(2), sz(3))/(zrange(2)-zrange(1))*sz(3);
end

function  RA = referencegen(sz,dp)
RA = imref3d(sz,[-dp(1,1)*sz(2)/2 dp(1,1)*sz(2)/2],...
    [-dp(1,2)*sz(1)/2 dp(1,2)*sz(1)/2],[dp(1,3) sz(3)*dp(1,3)]);
end

function reftform = pcalignrigid(filename,l1,l2,dp1,dp2,surfpoint1,...
    surfpoint2_ilm,surfpoint2_rpe)

%Find matching points
load(filename,'locs','matchedpoints') 
xi = locs(:,1); yi = locs(:,2);
x2i = matchedpoints(:,1); y2i = matchedpoints(:,2);


matchedheight1 =  interp2(l1,xi,yi);
matchedpoints1 = [(xi-size(l1,1)/2)*dp1(1),(yi-size(l1,2)/2)*dp1(2),matchedheight1*dp1(3)];
matchedheight2 =  interp2(l2,x2i,y2i);
matchedpoints2 = [(x2i-size(l1,1)/2)*dp2(1),(y2i-size(l1,2)/2)*dp2(2),matchedheight2*dp2(3)];

try
    tform2 = estgeotform3d(matchedpoints1,matchedpoints2,"rigid", ...
        "MaxDistance",800);
catch
    trans = matchedpoints2 - matchedpoints1;
    tform2 = rigidtform3d([0 0 0],nanmean(trans));
end



% Let's see how well these 2 surfaces merge
surfpoint2 = [surfpoint2_ilm;surfpoint2_rpe];

% Let's look at intersecting regions and take slab
k = convhull(surfpoint2(:,1),surfpoint2(:,2));
transsurf1 = transformPointsForward(tform2, surfpoint1);
in = inpolygon(transsurf1(:,1),transsurf1(:,2),surfpoint2(k,1),surfpoint2(k,2));
transsurf1 = transsurf1(in,:);

temp1 = transformPointsInverse(tform2,transsurf1);
[~,reftform] = genpartialcloud_2layers(temp1,tform2,surfpoint2_ilm,...
    surfpoint2_rpe);

end

function [refcloud,reftform] = genpartialcloud_2layers(subcloud,tform,aligncloud,aligncloud2)
% Transform to aligns reference frame
inliersPtCloud = transformPointsForward(tform, subcloud);
% Get heights at transposed points
F = scatteredInterpolant(aligncloud(:,1),aligncloud(:,2),aligncloud(:,3));
height =  F(inliersPtCloud(:,1),inliersPtCloud(:,2));
%  Transform to correct reference frame
refcloud = [inliersPtCloud(:,1),inliersPtCloud(:,2),height];
F = scatteredInterpolant(aligncloud2(:,1),aligncloud2(:,2),aligncloud2(:,3));
height =  F(inliersPtCloud(:,1),inliersPtCloud(:,2));
refcloud2 = [inliersPtCloud(:,1),inliersPtCloud(:,2),height];
refcloud = [refcloud;refcloud2];
reftform = pcregistericp(pointCloud(subcloud),pointCloud(refcloud),...
    'InitialTransform',tform);
end

function [box,boxedge] = pcbox(coords,dp)
%Point for every coordinate within original frame
% Initialize values
z = dp(3):dp(3):coords(3)*dp(3);
x = (-coords(1)-1)/2*dp(1):dp(1):(coords(1)+1)/2*dp(1);
y = (-coords(2)-1)/2*dp(2):dp(2):(coords(2)+1)/2*dp(2);
[X,Y,Z] = meshgrid(x,y,z);
box = [X(:),Y(:),Z(:)];
boxedge = [(-coords(1)-1)/2*dp(1), (-coords(2)-1)/2*dp(2), dp(3);
    (-coords(1)-1)/2*dp(1), (coords(2)+1)/2*dp(2), dp(3);
    (coords(1)+1)/2*dp(1), (-coords(2)-1)/2*dp(2), dp(3);
    (coords(1)+1)/2*dp(1), (coords(2)+1)/2*dp(2), dp(3);
    (-coords(1)-1)/2*dp(1), (-coords(2)-1)/2*dp(2), coords(3)*dp(3);
    (-coords(1)-1)/2*dp(1), (coords(2)+1)/2*dp(2), coords(3)*dp(3);
    (coords(1)+1)/2*dp(1), (-coords(2)-1)/2*dp(2), coords(3)*dp(3);
    (coords(1)+1)/2*dp(1), (coords(2)+1)/2*dp(2), coords(3)*dp(3)];
end

function box = pcboxgen(frame_sz,xcoord,ycoord,zcoord)
x = linspace(xcoord(1),xcoord(2),frame_sz(2));
y = linspace(ycoord(1),ycoord(2),frame_sz(1));
z = linspace(zcoord(1),zcoord(2),frame_sz(3));
[X,Y,Z] = meshgrid(x,y,z);
box = [X(:),Y(:),Z(:)];
end

function surfpoint = findsurf(height,dp) %Takes origin as center of FOV
idx = find(height>0);
[r,c] = ind2sub(size(height),idx);
r = r-(size(height,1)+1)/2; c = c-(size(height,2)+1)/2;
r = r*dp(2); c = c*dp(1);
surfpoint = [c,r,height(idx)*dp(3)];
end

function [refcloud,reftform] = genpartialcloud2(subcloud,tform,aligncloud_ilm, ...
    aligncloud_rpe)
% Transform to aligns reference frame
inliersPtCloud = transformPointsForward(tform, subcloud);
% Get heights at transposed points
F = scatteredInterpolant(aligncloud_ilm(:,1),aligncloud_ilm(:,2),aligncloud_ilm(:,3));
height =  F(inliersPtCloud(:,1),inliersPtCloud(:,2));
%  Transform to correct reference frame
refcloud = [inliersPtCloud(:,1),inliersPtCloud(:,2),height];
% Get heights at transposed points
F = scatteredInterpolant(aligncloud_rpe(:,1),aligncloud_rpe(:,2),aligncloud_rpe(:,3));
height =  F(inliersPtCloud(:,1),inliersPtCloud(:,2));
%  Transform to correct reference frame
temp = [inliersPtCloud(:,1),inliersPtCloud(:,2),height];
refcloud = [refcloud;temp];
[reftform, ~, rmse] = pcregistericp(pointCloud(subcloud),pointCloud(refcloud),...
    'InitialTransform',tform);
%Realign if the error in the reigstration is too great
realign = 0;
if rmse > 20
    realign = 1;
    iter = 1;
    reftform = pcregistericp(pointCloud(subcloud),pointCloud(refcloud),...
    'InitialTransform',tform,'InlierRatio',0.6);
end
end

function outImage = blendRegions(A,B)
outImage = cat(4,A,B);
outImage = max(outImage,[],4);
end


