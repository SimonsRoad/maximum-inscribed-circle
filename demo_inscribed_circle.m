clc;clear;close all

%I=imread('hand_contour.png');
%I = imread('weird_region.png');
I = imread('11.bmp');

nhood = true(9);
S = stdfilt(I,nhood);
gray = mat2gray(S);
if length(size(gray)) == 2
    umb = gray;
else
    umb = rgb2gray(gray);
end
level = graythresh(umb);
bw = imbinarize(umb,level);

% se = strel('disk',5);
% se1 = strel('disk',2);
% bw0 = imerode(imdilate(bw,se),se1);

BW = imfill(bw,'holes');
result = bwareaopen(BW,800);

[boundaries,L] = bwboundaries(result);

if length(boundaries) < 1
    disp('no boundary detected...')
    return;
else % 保留提取到的最大轮廓
    blen = zeros(1,length(boundaries));
    for i = 1:length(boundaries)
        blen(i) = size(boundaries{i},1);
    end
    [~,cidx] = max(blen);
    boundary=boundaries{cidx};
end

% 边界等于1，其余为0
contourImage = logical(zeros(size(result,1),size(result,2))==1);
for i = 1:size(boundary,1)
    rowXY = boundary(i,:);
    contourImage(rowXY(1),rowXY(2)) = 1;
end

%% 方法一
[r1,cx1,cy1] = inscribed_circle_poly(contourImage,1);


%% 方法二，与一类似
%[r2,cx2,cy2] = inscribed_circle_poly2(contourImage,1);

%% 方法三
%[r3,cx3,cy3] = inscribed_circle_skel(result,1);

%% 方法四
[r4,cx4,cy4] = inscribed_circle_opti(contourImage,1);

% inpoly test
% AP = [10,10;20,30;30,20;40,10;30,0;20,10];
% TP = [20,20;20,30];
% figure
% plot(AP(:,1),AP(:,2),'k-')
% hold on
% scatter(TP(:,1),TP(:,2),'filled')
% iores = inpoly([TP(:,1),TP(:,2)],[AP(:,1),AP(:,2)])
disp('bingo...')