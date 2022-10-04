function [cx,cy,r]=inscribed_circle_poly2(contourImage,display)
% 输入为原始图像轮廓图
% given a polygon [x,y] find the inner circle [cx,cy,r]
% by Tolga Birdal
% A sample application and a function for solving the maximum inscribed
% circle problem. 
% Unlike my other submission "Maximum Inscribed Circle using Distance 
% Transform", this algorithm is subpixel accurate. It operates only on the
% polygon and not the image points. Therefore, if the polygon is given in
% sub-pixels, the result will be accurate. 
% I use an O(n log(n)) algorithm as follows:
% Construct the Voronoi Diagram of the polygon.
% For Voronoi nodes which are inside the polygon:
%       Find the node with the maximum distance to edges in P. This node is
%       the centre of the maximum inscribed circle.
% 
% For more details on the problem itself please checkout my previous 
% submission as mentioned above.
% 
% To speed things up:
% Replace "inpolygon" function by Bruno Lunog's faster implementation:
% "2D polygon interior detection" :
% http://www.mathworks.com/matlabcentral/fileexchange/27840-2d-polygon-inte
% rior-detection
% Copyright (c) 2011, Tolga Birdal <http://www.tbirdal.me>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
[y,x]=find(contourImage==1,1, 'first');
contour = bwtraceboundary(contourImage, [y(1),x(1)], 'E', 8);
x=contour(:,2);
y=contour(:,1);

% make a voronoi diagram
[vx,vy]=voronoi(x,y);

% find voronoi nodes inside the polygon [x,y]
Vx=vx(:);
Vy=vy(:);
figure
imshow(contourImage)
hold on
plot(Vx,Vy,'+')

% Here, you could also use a faster version of inpolygon
%IN=inpolygon(Vx,Vy, x,y); % 轮廓内的所有点
IN=inpoly([Vx,Vy], [x,y]); % 轮廓内的所有点
ind=find(IN==1);
Vx=Vx(ind);
Vy=Vy(ind);

% maximize the distance of each voronoi node to the closest node on the
% polygon.
minDist=0;
minDistInd=-1;
for i=1:length(Vx)
    dx=(Vx(i)-x);
    dy=(Vy(i)-y);
    r=min(dx.*dx+dy.*dy);
    if (r>minDist)
        minDist=r;
        minDistInd=i;
    end
end

% take the center and radius
cx=Vx(minDistInd);
cy=Vy(minDistInd);
r=sqrt(minDist);
toc

% display result
if (~exist('display','var'))
    display=1;
end
if (display)
    figure
    imshow(contourImage,[]);
    hold on
    plot(x,y,'m','LineWidth',0.2);
    
    plot(cx,cy,'g-*')
    theta = [linspace(0,2*pi) 0];
    hold on
    plot(cos(theta)*r+cx,sin(theta)*r+cy,'color','r','LineWidth', 0.1);
end

end
