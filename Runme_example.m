% Base home directory with all the data
homedir = 'E:\RetinalSLAM\ExampleDemo_Copy\';

% Image number of each data acquistion. Folder 632 contains the b-scans of
% volume 1 and their respective segmentations. 634 contains the data for
% volume 2 etc.
framenum = [632, 634, 636, 638, 640, 644, 646, 648, 650];

% Overlapping volumes denoted in startcoord and endcoord. In this example,
% volume 2 (folder 634) overlaps with volume 1 (folder 632). Thus,
% startcoord(1) = 2 and endcoord(1) = 1. 
startcoord = [2,3,4,5,6,7,8,9];
endcoord = [1,1,1,1,5,2,2,2];

% Assigned file name to the matching points. R634_632 saves the matching
% point between volume 2 (folder 634) and volume 1 (folder 632). The name
% of these files is not important so long as they are different
coordfiles = ["R634_632.mat";"R636_632.mat";"R638_632.mat";
    "R640_632.mat";"R644_640.mat";"R646_634.mat";
    "R648_634.mat";"R650_634.mat"];

% Set to false if you don't want to do derformable registration.
deformable = true;

% Coordinate system used for volumetric montaging. In this case, we are
% using the coordinate system of volume 1 (folder 632). Recommend choosing
% volume with the most overlapping volumes although this can easily be
% altered
refvol = 1;

% Pixel size of x (6000/360 microns), y (6000/360 microns), and 
% z(2.6 microns) in the original acquistion 
dp = repmat([6000/360 6000/360 2.6],length(framenum),1);

% Pixel size of final reconstruction
sx = [30, 30, 5];

% Do the point matching (Step 3 in paper) between en-faces of volumes
savecond = true;
for ii = 1:length(startcoord)
    % First volume 
    load([homedir num2str(framenum(startcoord(ii))) '\frame_3D.mat'],'frame_3D')
    frame1 = frame_3D;
    % 2nd volume to match with 1st volume
    load([homedir num2str(framenum(endcoord(ii))) '\frame_3D.mat'],'frame_3D')
    frame2 = frame_3D;
    % Saves matching points between 1st and 2nd volumes in filename with
    % 3rd input
    automatedmatching_quadrants(frame1,frame2,char(coordfiles(ii)),savecond)
end

% Do the actual volumetric montaging
F = volume_merge3D(homedir,framenum,coordfiles,startcoord,endcoord,...
    deformable,refvol,dp,sx);

volumeSegmenter(F)
enface = squeeze(mean(maxk(F,15,3),3));
imagesc(enface)