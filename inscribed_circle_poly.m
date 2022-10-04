function [R,cx,cy]=inscribed_circle_poly(ContourImage, display)
% --- 函数说明 ---
% 获取任意图形的最大内切圆
%
% --- 输入 ---
% ContourImage : 目标图像轮廓图
% display : 是否输出图形结果展示
%
% --- 输出 ---
% R : 最大内切圆半径 
% cx : 内切圆圆心x坐标 
% cy : 内切圆圆心y坐标 
%
%     by Tolga Birdal
% 
% 
%     Maximum Inscribed Circle
% 
%     Or in other words, "largest inner circle" , "maximum empty circle" etc.
% 
%     This is a very common problem in computational geometry, and it is not 
%     simple to solve efficiently.
%     Addressing 2D image/contour processing, I couldn't find a good 
%     implementation on the web. Generally, the reasonable way of solving
%     this problem is to make use of Voronoi Diagrams, which are generally 
%     O(nlogn). 
% 
%     After analyzing the problem a bit, I noticed that the problem can 
%     easily and approximately be solved using well-known distance transform. 
% 
%     Here is how:
% 
%     The computational aim can be written as:  
%     (x, y) maximizes r = min_{i} r_{i} 
%     where r_i = ||(x_i, y_i) ? (x, y)|| and d_i = r_i ? r
%     
%     (x_i, y_i): Pairs data points
%     (x, y), r : Pair, scalar circle centre and radius
%     
%     In non-mathematical terms: 
%     
%     1.  The center of the maximum inscribed circle will lie inside the 
%         polygon
%     2.  The center of such a circle will be furthest from any point on the 
%         edges of the polygon.
% 
%     So we seek for the point that lies inside the polygon and has maximal 
%     distance to the closest edge. This is exactly the maximum pixel of the
%     distance transform lying inside the contour.
% 
%     Notice that this approach is completely in-precise. It is only 
%     pixel-precise and never subpixel accurate. However, unlike 
%     optimization approaches, it does guarantee a global convergence. In 
%     the case of ambiguity, any of the solutions will be valid. 
% 
%     To detect the points inside the region, inpolygon remains very slow. 
%     So, I make use of the great code of Darren Engwirda, here. As well as 
%     being contained in this package, it can also be downloaded from:
%     http://www.mathworks.com/matlabcentral/fileexchange/10391-fast-points-in-polygon-test
% 
%     Here are other implemnatations, which are more accurate, but much slower
%     than my approach (only slower in Matalb of course!)
%     
%     Using 
%     http://www.mathworks.com/matlabcentral/fileexchange/2794-cvoronoi
% 
%     Using "Largest pixel that doesn't cross any point" approach:
%     http://www.mathworks.com/matlabcentral/newsreader/view_thread/283296
%   
%     Here is a sample call: 
%
%       I=imread('hand_contour.png');
%       [R cx cy]=max_inscribed_circle(I)       
%
%   Cheers,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the contour
% figure
% imshow(ContourImage)

tic;
sz=size(ContourImage);
[py,px]=find(ContourImage==1,1, 'first'); %边界等于1
% hold on
% plot(px,py,'b*')
contour = bwtraceboundary(ContourImage, [py(1), px(1)], 'W', 8); % 搜索轮廓的所有点
X=contour(:,2);
Y=contour(:,1);

% find the maximum inscribed circle:
% The point that has the maximum distance inside the given contour is the
% center. The distance of to the closest edge (tangent edge) is the radius.

BW=bwdist(logical(ContourImage));
[Mx, My]=meshgrid(1:sz(2), 1:sz(1));
% figure
% imshow(ContourImage)
% hold on
% plot(Mx,My,'+')

[Vin,Von]=inpoly([Mx(:),My(:)],[X,Y]);
%[Vin,Von]=inpolygon(Mx(:),My(:),X,Y);
ind=sub2ind(sz, My(Vin),Mx(Vin));
[R,RInd]=max(BW(ind)); 
R=R(1); RInd=RInd(1); % handle multiple solutions: Just take first.
[cy,cx]=ind2sub(sz, ind(RInd));
toc

% display result
if (~exist('display','var'))
    display=1;
end
if (display)
    BW(cy,cx)=R+20; % to emphasize the center
    figure,imshow(BW,[]);
    hold on
    plot(X,Y,'m','LineWidth',0.2);
    
    plot(cx,cy,'g-*');
    theta = [linspace(0,2*pi) 0];
    hold on
    plot(cos(theta)*R+cx,sin(theta)*R+cy,'color','r','LineWidth', 0.1);
end

end
