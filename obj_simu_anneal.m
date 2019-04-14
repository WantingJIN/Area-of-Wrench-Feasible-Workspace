function workspace = obj_simu_anneal( x )

workspace = 0 ;
%x0 contains the set of variables: [b,L,l0,k1,k2,Fmax]


b = 1;
L = b*x(1);
l0=L-b;
% k1 = x(2);
% k2 = x(3);
K = 20*L;
k1=1/(1+x(2))*K ;
k2=x(2)/(1+x(2))*K;

Fmax = 100;
%joint limits
lb1=-pi/2;
ub1=pi/2;
lb2=-pi/2;
ub2=pi/2;



k11=k1;
k12=k1;
k21=k2;
k22=k2;

%% compute (x,y) for all values of (theta1,theta2) : direct model
%generating theta 1 and theta2
n = 200;
TH1=linspace(lb1,ub1,n);
TH2=linspace(lb2,ub2,n);

% X=zeros(length(TH1),length(TH2));
% Y=zeros(length(TH1),length(TH2));
X=[];
Y=[];
TH=[];
for i=1:length(TH1)
    for j=1:length(TH2)
    theta1=TH1(i);
    theta2=TH2(j);
    
    TH=[TH;theta1,theta2];
    
    end
end
P=zeros(size(TH));
for i=1:size(TH,1)
    
    theta1=TH(i,1);
    theta2=TH(i,2);
    %reference point (x,y)
    P(i,:)=[-sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)-sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2+2*theta1),sqrt(L^2-b^2*cos(theta1)^2)*cos(theta1)+sqrt(L^2-b^2*cos(theta2)^2)*cos(theta2+2*theta1)];

end

%% compute F1,F2 for all values of theta1, theta2 

F=zeros(size(TH));
Work=zeros(size(TH));
Fmin=0;
Fmax=100;

for i=1:size(TH,1)
    
    theta1=TH(i,1);
    theta2=TH(i,2);
    
    F1 = (cos(theta2)^2*sin(theta1)*b^3*k21-cos(theta2)^2*sin(theta1)*b^3*k22-sin(theta2)*cos(theta1)^2*b^3*k11+sin(theta2)*cos(theta1)^2*b^3*k12+.5*L^2*sin(theta2)*b*k11-.5*L^2*sin(theta2)*b*k12-.5*L^2*sin(theta1)*b*k21+.5*L^2*sin(theta1)*b*k22-.5*sin(theta1)*b^3*k21+.5*sin(theta1)*b^3*k22+.5*sin(theta2)*b^3*k11-.5*sin(theta2)*b^3*k12+.5*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k22+.5*sqrt(L^2-b^2*cos(theta2)^2)*b^2*k11-.5*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k21-.5*sqrt(L^2-b^2*cos(theta2)^2)*b^2*k12-.5*L^2*sqrt(L^2-b^2*cos(theta2)^2)*k12+.5*L^2*sqrt(L^2-b^2*cos(theta2)^2)*k11-.5*L^2*sqrt(L^2-b^2*cos(theta1)^2)*k21+.5*L^2*sqrt(L^2-b^2*cos(theta1)^2)*k22-.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k22*l0+sqrt(L^2-b^2*cos(theta2)^2)*cos(theta1)^2*b^2*k12+cos(theta2)^2*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k21-.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k22*l0-1.*sqrt(L^2-b^2*cos(theta2)^2)*cos(theta1)^2*b^2*k11+.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k21*l0-1.*cos(theta2)^2*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k22-.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k11*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k12*l0-1.*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b^2*k12+sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k21+.5*sin(theta2)*sin(theta1)*b^2*k11*l0+.5*sin(theta2)*sin(theta1)*b^2*k12*l0-.5*sin(theta2)*sin(theta1)*b^2*k21*l0-.5*sin(theta2)*sin(theta1)*b^2*k22*l0+sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k22-1.*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b^2*k11+sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sin(theta1)*b^2*k21+sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sin(theta1)*b^2*k22-.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k21*l0+.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k12*l0-.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k11*l0-.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k22*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k12*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k21*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k11*l0-1.*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b*k12-1.*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b*k11)/((sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)-sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2))*b);
    F2 = (-cos(theta2)^2*sin(theta1)*b^3*k21+cos(theta2)^2*sin(theta1)*b^3*k22+sin(theta2)*cos(theta1)^2*b^3*k11-sin(theta2)*cos(theta1)^2*b^3*k12-.5*L^2*sin(theta2)*b*k11+.5*L^2*sin(theta2)*b*k12+.5*L^2*sin(theta1)*b*k21-.5*L^2*sin(theta1)*b*k22+.5*sin(theta1)*b^3*k21-.5*sin(theta1)*b^3*k22-.5*sin(theta2)*b^3*k11+.5*sin(theta2)*b^3*k12+.5*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k22+.5*sqrt(L^2-b^2*cos(theta2)^2)*b^2*k11-.5*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k21-.5*sqrt(L^2-b^2*cos(theta2)^2)*b^2*k12-.5*L^2*sqrt(L^2-b^2*cos(theta2)^2)*k12+.5*L^2*sqrt(L^2-b^2*cos(theta2)^2)*k11-.5*L^2*sqrt(L^2-b^2*cos(theta1)^2)*k21+.5*L^2*sqrt(L^2-b^2*cos(theta1)^2)*k22-.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k22*l0+sqrt(L^2-b^2*cos(theta2)^2)*cos(theta1)^2*b^2*k12+cos(theta2)^2*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k21-.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k22*l0-1.*sqrt(L^2-b^2*cos(theta2)^2)*cos(theta1)^2*b^2*k11+.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k21*l0-1.*cos(theta2)^2*sqrt(L^2-b^2*cos(theta1)^2)*b^2*k22-.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k11*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*k12*l0+sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b^2*k12+sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k21-.5*sin(theta2)*sin(theta1)*b^2*k11*l0-.5*sin(theta2)*sin(theta1)*b^2*k12*l0+.5*sin(theta2)*sin(theta1)*b^2*k21*l0+.5*sin(theta2)*sin(theta1)*b^2*k22*l0+sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k22+sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b^2*k11-1.*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sin(theta1)*b^2*k21-1.*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta2)*sin(theta1)*b^2*k22-.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k21*l0-.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k12*l0+.5*sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2)*b*k11*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k22*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k12*l0-.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k21*l0+.5*sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)*b*k11*l0-1.*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b*k12-1.*sqrt(L^2-b^2*cos(theta2)^2)*sqrt(L^2-b^2*cos(theta1)^2)*sin(theta1)*b*k11)/((sqrt(L^2-b^2*cos(theta2)^2)*sin(theta1)-sin(theta2)*sqrt(L^2-b^2*cos(theta1)^2))*b);
     F(i,:)=[F1,F2];
    if F(i,1)>Fmin && F(i,2)>Fmin && F(i,1)<Fmax && F(i,2)<Fmax 
        Work(i) = 1;
    else
        Work(i) = NaN;

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
% N=100;
% M=100;
% xdot=(xmax-xmin)/N;
% ydot=(ymax-ymin)/M;
xdot=0.01;
ydot=0.01;
N=round((xmax-xmin)/xdot);
M=round((ymax-ymin)/ydot);
map=zeros(N+1,M+1);
list=[];
 if (xdot>0) && (ydot>0)
    for i=1:size(TH,1)
        if Work(i)==1
             list=[list i];
            map_x = fix((P(i,1)-xmin)/xdot) + 1;
            map_y = fix((P(i,2)-ymin)/ydot) + 1;
            map(map_x,map_y) = map(map_x,map_y) + 1;
        end
    end
    countspace=0;
    for i=1:N
        for j=1:M
            if map(i,j)>0
                countspace=countspace+1;
            end
        end
    end
    Space=xdot*ydot*countspace*10000;
 else Space=0;
 end
% workspace = -countspace/size(list,2)*100
workspace = -Space
end


