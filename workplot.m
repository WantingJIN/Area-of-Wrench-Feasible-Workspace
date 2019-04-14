%% Initialisation
clear all
close all
tic
%global x
% x=[4.8166 18.5436 0.2976 0.0115 0.3285];
% x=[1 3 0.5 20 30];
x= [1 20 10 0.2 0.3];
b = x(1);
L = x(2);
l0 = x(3);


%joint limits
lb1=-pi/2;
ub1=pi/2;
lb2=-pi/2;
ub2=pi/2;

%spring stiffness
k1 = x(4);
k11=k1;
k12=k1;
k2 = x(5);
k21=k2;
k22=k2;

%% compute (x,y) for all values of (theta1,theta2) : direct model
%generating theta 1 and theta2
% TH1=linspace(lb1,ub1,n);
% TH2=linspace(lb2,ub2,n);
% 
% % X=zeros(length(TH1),length(TH2));
% % Y=zeros(length(TH1),length(TH2));
% X=[];
% Y=[];
% TH=[];
% 
% for i=1:length(TH1)
%     for j=1:length(TH2)
%     theta1=TH1(i);
%     theta2=TH2(j);
%    
%     TH=[TH;theta1,theta2];
%     
%     end
% end
N=200000;
n=100000;
rng(200);
TH_wide=-pi/2 + pi.*rand(N,2);
TH_narrow=-pi/14 + pi/7.*rand(n,2);
TH=[TH_wide;TH_narrow];
P=zeros(size(TH));
% P is the x,y coordinate of top bar center
for i=1:size(TH,1)
    
    theta1=TH(i,1);
    theta2=TH(i,2);
    x=-(1/2)*sqrt(-2.*b^2*cos(2*theta1)+4*L^2-2*b^2)*sin(theta1)-(1/2)*sqrt(-2.*b^2*cos(2*theta2)+4*L^2-2*b^2)*sin(theta2+2*theta1);
    y=(1/2)*sqrt(-2.*b^2*cos(2*theta1)+4*L^2-2*b^2)*cos(theta1)+(1/2)*sqrt(-2.*b^2*cos(2*theta2)+4*L^2-2*b^2)*cos(theta2+2*theta1);
%     P(i,:)=[-sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)-sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2+2*theta1),sqrt(L^2-b^2*cos(theta1)^2)*cos(theta1)+sqrt(L^2-b^2*cos(theta2)^2)*cos(theta2+2*theta1)];
    P(i,:)=[x,y];
end

%% compute F1,F2 for all values of theta1, theta2 

F=zeros(size(TH));
Fmin=0;
Fmax=100;

for i=1:size(TH,1)
    
    theta1=TH(i,1);
    theta2=TH(i,2);

    F1 = (cos(theta2)^2*sin(theta1)*b^3*k21-cos(theta2)^2*sin(theta1)*b^3*k22-sin(theta2)*cos(theta1)^2*b^3*k11+sin(theta2)*cos(theta1)^2*b^3*k12+.5*L^2*sin(theta2)*b*k11-.5*L^2*sin(theta2)*b*k12-.5*L^2*sin(theta1)*b*k21+.5*L^2*sin(theta1)*b*k22-.5*sin(theta1)*b^3*k21+.5*sin(theta1)*b^3*k22+.5*sin(theta2)*b^3*k11-.5*sin(theta2)*b^3*k12+.5*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k22+.5*sqrt(L^2-b^2*cos(theta2)^2)*b^2*k11-.5*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k21-.5*sqrt(L^2-b^2*cos(theta2)^2)*b^2*k12-.5*L^2*sqrt(L^2-b^2*cos(theta2)^2)*k12+.5*L^2*sqrt(L^2-b^2*cos(theta2)^2)*k11-.5*L^2*sqrt(L^2-b^2*cos(theta1)^2)*k21+.5*L^2*sqrt(L^2-b^2*cos(theta1)^2)*k22-.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k22*l0+sqrt(L^2-b^2*cos(theta2)^2)*cos(theta1)^2*b^2*k12+cos(theta2)^2*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k21-.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k22*l0-1.*sqrt(L^2-b^2*cos(theta2)^2)*cos(theta1)^2*b^2*k11+.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k21*l0-1.*cos(theta2)^2*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k22-.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k11*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k12*l0-1.*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b^2*k12+sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k21+.5*sin(theta2)*sin(theta1)*b^2*k11*l0+.5*sin(theta2)*sin(theta1)*b^2*k12*l0-.5*sin(theta2)*sin(theta1)*b^2*k21*l0-.5*sin(theta2)*sin(theta1)*b^2*k22*l0+sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k22-1.*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b^2*k11+sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sin(theta1)*b^2*k21+sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sin(theta1)*b^2*k22-.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k21*l0+.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k12*l0-.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k11*l0-.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k22*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k12*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k21*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k11*l0-1.*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b*k12-1.*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b*k11)/((sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)-sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2))*b);
    F2 = (-cos(theta2)^2*sin(theta1)*b^3*k21+cos(theta2)^2*sin(theta1)*b^3*k22+sin(theta2)*cos(theta1)^2*b^3*k11-sin(theta2)*cos(theta1)^2*b^3*k12-.5*L^2*sin(theta2)*b*k11+.5*L^2*sin(theta2)*b*k12+.5*L^2*sin(theta1)*b*k21-.5*L^2*sin(theta1)*b*k22+.5*sin(theta1)*b^3*k21-.5*sin(theta1)*b^3*k22-.5*sin(theta2)*b^3*k11+.5*sin(theta2)*b^3*k12+.5*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k22+.5*sqrt(L^2-b^2*cos(theta2)^2)*b^2*k11-.5*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k21-.5*sqrt(L^2-b^2*cos(theta2)^2)*b^2*k12-.5*L^2*sqrt(L^2-b^2*cos(theta2)^2)*k12+.5*L^2*sqrt(L^2-b^2*cos(theta2)^2)*k11-.5*L^2*sqrt(L^2-b^2*cos(theta1)^2)*k21+.5*L^2*sqrt(L^2-b^2*cos(theta1)^2)*k22-.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k22*l0+sqrt(L^2-b^2*cos(theta2)^2)*cos(theta1)^2*b^2*k12+cos(theta2)^2*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k21-.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k22*l0-1.*sqrt(L^2-b^2*cos(theta2)^2)*cos(theta1)^2*b^2*k11+.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k21*l0-1.*cos(theta2)^2*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k22-.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k11*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k12*l0+sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b^2*k12+sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k21-.5*sin(theta2)*sin(theta1)*b^2*k11*l0-.5*sin(theta2)*sin(theta1)*b^2*k12*l0+.5*sin(theta2)*sin(theta1)*b^2*k21*l0+.5*sin(theta2)*sin(theta1)*b^2*k22*l0+sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k22+sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b^2*k11-1.*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sin(theta1)*b^2*k21-1.*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sin(theta1)*b^2*k22-.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k21*l0-.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k12*l0+.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k11*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k22*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k12*l0-.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k21*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k11*l0-1.*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b*k12-1.*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b*k11)/((sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)-sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2))*b);
    F(i,:)=[F1,F2];

end


%% Stiffness computation

count = 0;
for i=1:size(TH,1)
    
    theta1=TH(i,1);
    theta2=TH(i,2); 
    F1=F(i,1);
    F2=F(i,2);
    %Compute the possible position if 0<F<Fmax
if F(i,1)>Fmin && F(i,2)>Fmin && F(i,1)<Fmax && F(i,2)<Fmax    
    count=count+1;
    Work(i)=1;
else
    Work(i)= NaN;
end

end
xmin=10;
xmax=-10;
ymin=10;
ymax=-10;

for i=1:size(TH,1)
    if Work(i)==1
       if P(i,1) > xmax
           xmax = P(i,1);
       end
       if P(i,1) < xmin
           xmin = P(i,1);
       end
       if P(i,2) > ymax
           ymax = P(i,2);
       end
       if P(i,2) < ymin
           ymin = P(i,2);
       end
    end
end
% xmax=max(P(:,1));
% xmin=min(P(:,1));
% ymax=max(P(:,2));
% ymin=min(P(:,2));
N=500;
M=500;
xdot=(xmax-xmin)/N;
ydot=(ymax-ymin)/M;
map=zeros(N+1,M+1);
% xdot=0.05;
% ydot=0.05;
% N=round((xmax-xmin)/xdot);
% M=round((ymax-ymin)/ydot);
% map=zeros(N+1,M+1);
list=[];
for i=1:size(TH,1)
    if Work(i)==1
        list=[list i];
        map_x = fix((P(i,1)-xmin)/xdot) + 1;
        map_y = fix((P(i,2)-ymin)/ydot) + 1;
        map(map_x,map_y) = map(map_x,map_y) + 1;
    end
end

%find edge
% 
edge=zeros(fix((N+1)/2)+1,M+1);
area=zeros(fix((N+1)/2)+1,M+1);
edgex=zeros(fix((N+1)/2)+1,M+1);
edgey=zeros(fix((N+1)/2)+1,M+1);
map=map(1:fix((N+1)/2),:);
for i= 1:fix((N+1)/2) %% the right side and the left side are symmetic
        edger_x=i;       %find the right edge
        edger_y=max(find(map(i,:)));  
        edgey(edger_x,edger_y)=1;
        edgel_x=i;       %find the left edge
        edgel_y=min(find(map(i,:)));  
        edgey(edgel_x,edgel_y)=1;
%         imshow(rot90(edge,1))
        if(edger_y-edgel_y>8*xdot)
             middle_y=fix((edger_y-edgel_y)/2)+edgel_y;
             if map(i,middle_y)==0    % If the middle point doesn't in workspace, then there is a hole in it
                edger_x2=i;
                edger_y2=min(find(map(i,middle_y:M)));
                edgey(edger_x2,edger_y2+middle_y-1)=1;
                edgel_x2=i;       %find the left edge
                edgel_y2=max(find(map(i,1:middle_y)));  
                edgey(edgel_x2,edgel_y2)=1;
             end

        end
%         imshow(rot90(edgey,1))
end
figure(1)
imshow(rot90(edgey,1))
title('The vertical edge')
for j=1:M
        edgeu_y=j;       %find the up edge
        edgeu_x=max(find(map(:,j)));  
        edgex(edgeu_x,edgeu_y)=1;
        edgeb_y=j;       %find the bottom edge
        edgeb_x=min(find(map(:,j)));  
        edgex(edgeb_x,edgeb_y)=1;
%         imshow(rot90(edgex,1))
        if(edgeu_x-edgeb_x>100*ydot)
           middle_x=fix((edgeu_x-edgeb_x)/2)+edgeb_x;
           if map(middle_x,j)==0
                edgeu_x2=max(find(map(1:middle_x,j)));
                edgeu_y2=j;
                edgex(edgeu_x2,edgeu_y2)=1;
                edgeb_x2=min(find(map(middle_x:fix((N+1)/2),j)));       
                edgeb_y2=j; 
                edgex(edgeb_x2+middle_x,edgeb_y2)=1;
%                 imshow(rot90(edgex,1))

           end
        end
end
figure(2)
imshow(rot90(edgex,1))
title('The lateral edge')
figure(3)
edge=edgex+edgey;
imshow(rot90(edge,1))
title('The whole edge')
% edge=edgey;

for i=1:fix((N+1)/2)
    edy=find(edgey(i,:));
    if size(edy,2)==4
        if edy(2)-edy(1)<M/2
            edgey(i,edy(1):edy(2))=1;
        else
            if map(i,fix((edy(2)-edy(1))/2))>0
                edgey(i,edy(1):edy(2))=1;
            end
        end
        if edy(4)-edy(3)<M/2
            edgey(i,edy(3):edy(4))=1;
        else
            if map(i,fix((edy(4)-edy(3))/2))>0
                edgey(i,edy(3):edy(4))=1;
            end
        end
    end
    if size(edy,2)==3||size(edy,2)==2
        if edy(2)-edy(1)<M/2
        edgey(i,edy(1):edy(2))=1;
        else
            if map(i,fix((edy(2)-edy(1))/2))>0
                edgey(i,edy(1):edy(2))=1;
            end
        end
    end
%     imshow(rot90(edgey,1))
end
figure(4)
imshow(rot90(edgey,1))
title('The WFW fill by the vertical edge')
for j=1:M
    edx=find(edgex(:,j));
    if size(edx,1)==4
        if edx(2)-edx(1)<M/2
        edgex(edx(1):edx(2),j)=1;
        else
            if map(fix((edx(2)-edx(1))/2),j)>0
                edgex(edx(1):edx(2),j)=1;
            end
            
        end
        if edx(4)-edx(3)<M/2
            edgex(edx(3):edx(4),j)=1;
        else
            if map(fix((edx(4)-edx(3))/2),j)>0
                edgex(edx(3):edx(4),j)=1;
            end
        end
    end
    
    if size(edx,1)==3||size(edx,1)==2
        if edx(2)-edx(1)<M/2
        edgex(edx(1):edx(2),j)=1;
        else
            if map(fix((edx(2)-edx(1))/2),j)>0
                edgex(edx(1):edx(2),j)=1;
            end
        end
    end
%     imshow(rot90(edgex,1))
end
figure(5)
imshow(rot90(edgex,1))
title('The WFW fill by the lateral edge')

figure(6)
subplot(1,3,1)
imshow(rot90(map,1))
title('Using adpative mesh')
subplot(1,3,2)
imshow(rot90(edge,1))
title('The edge of the workspace ')
area=edgex+edgey;
subplot(1,3,3)
imshow(rot90(area,1))
title('The whole area of the workspace ')

squarenumber=0;
for i=1:fix((N+1)/2)
    for j=1:M
        if area(i,j)>0
            squarenumber=squarenumber+1;
        end
    end
end
Space=xdot*ydot*squarenumber
%Plot 
% figure(2)
% sz=10;
% 
% subplot(1,2,1)
% hold on
% axis equal
% scatter(TH(:,1),TH(:,2),sz,Work,'filled');
% subplot(1,2,2)
% hold on
% axis equal
% scatter(P(:,1),P(:,2),sz,Work,'filled');
% %title('The wrench feasible workspce')
% toc