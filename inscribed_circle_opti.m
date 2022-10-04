function [r,cx,cy] = inscribed_circle_opti(contourImage,display)
% --- 函数说明 ---
% 获取任意图形的最大内切圆
% 求解非线性规划问题，目标函数为：max r(即半径)
% 约束条件：min(sqrt((x-a)^2+(y-b)^2))>=r
% 其中(x,y)为图形边界上点的坐标，(a,b)为圆心坐标，r为圆半径。
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
% Programmer: Robin An, 2021-10-01
% last modified by Robin An on 2021-10-01
% reference: http://www.matlabsky.com/thread-39620-1-1.html
%            天津科技大学数学系 谢中华
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fun = @(x)-x(3); %x = beta = [beta(1),beta(2),beta(3)]
%BW = imread(filename);
[y,x] = find(contourImage);
warning off

% fun 目标函数
% funcontest约束条件
% fmincon(FUN,X,A,B,Aeq,Beq,LB,UB,NONLCON,options,varargin)
initVal = [mean(x),mean(y),50];
lb = [min(x);min(y);0];
ub = [max(x);max(y);min(max(y)-min(y),max(x)-min(x))];%横纵坐标差值最大的最小一个
[beta,fval] = fmincon(fun,initVal,[],[],[],[],lb,ub,@funcontest);
% fun的输入为beta（beta(1)、beta(2)、beta(3)）

cx = beta(1);
cy = beta(2);
r = beta(3);
toc

% res = inpoly([cx,cy],[x,y])
% res1 = inpoly_se(cx,cy,x,y)
% res2 = inpoly_bl(cx,cy,x,y)
% res3 = inpoly([406,384],[x,y])
% res4 = inpoly_se(406,384,x,y)
% res5 = inpoly_bl(406,384,x,y)
% figure
% imshow(contourImage)
% hold on
% plot(cx,cy,'g-*') % 282,236
% hold on
% plot(406,384,'b-+')

% display result
if (~exist('display','var'))
    display=1;
end
if (display)
    figure
    imshow(contourImage)
    hold on
    plot(beta(1),beta(2),'g-*')
    hold on
    
    t = linspace(0,2*pi,60);
    xd = beta(1)+beta(3)*cos(t);
    yd = beta(2)+beta(3)*sin(t);
    plot(xd,yd,'r')
    hold off
end

% 非线性约束子函数，C最大化
    function [C,Ceq] = funcontest(beta)
        % a = beta(1),b = beta(2),r = beta(3)
                
        C = beta(3)-min(sqrt((beta(1)-x).^2+(beta(2)-y).^2));
        
        % 点在轮廓内部
        %Ceq = [];
        Ceq = inpoly([beta(1),beta(2)],[x,y]) - 1;
    end
end
