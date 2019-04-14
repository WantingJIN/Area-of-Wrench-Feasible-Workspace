function workspace = scan_obj_new(x)
b = x(1);
L = x(2);
l0 = x(3);
k1 = x(4);
k2 = x(5);
Fmax = 100;

k11=k1;
k12=k1;
k21=k2;
k22=k2;

%define the area for the possible workspace
xmax=2*L;
xmin=-2*L;
ymax=2*sqrt(L^2-b^2);
ymin=-2*sqrt(L^2-b^2);

xdot=0.1;
ydot=0.1;
N=round((xmax-xmin)/xdot);
M=round((ymax-ymin)/ydot);
Space=0;
for i=1:N
    for j=1:M
        work=0;
        x=xmin+i*xdot;
        y=ymin+i*ydot;
        work=solveIGM(L,b,x,y);
        if(work)
            Space=Space+1;
        end
    end
end

% workspace = -countspace/size(list,2)*100
workspace = -Space
end