
T = 10;% init temperature
Tf = 0.01;% Final value
eps_ratioLb = 0.5;  % Tune eps big in order to find the global minimum 
% eps_L = 1;
% eps_ratio = 1.5;  %L/b
% eps_ratiol0b = 0.1;
% eps_k1 = 0.1;
eps_k2k1 = 0.5;
%eps_F = 5;


% eps = [eps_b;eps_L;eps_l0;eps_k1;eps_k2];
eps = [eps_ratioLb;eps_k2k1];
maxiter = 10;
x =  [1.1;4]; %guess value
X_opt = x;
f_opt = obj_simu_anneal(X_opt);

while T > Tf
    iter = 0;
    while iter < maxiter
        y = 2.*eps.*rand(2,1) + x - eps;
%         if y(1)<0||y(2)<0||y(3)<0||y(4)<0||y(1)<1||y(1)>3
        if y(1)<1||y(2)<0
            break;
        end
        Fy = obj_simu_anneal(y);
        Fx = obj_simu_anneal(x);
        delta_E = Fy - Fx;
        if delta_E < 0
            x = y;
            if Fx < f_opt
                X_opt = y;
                f_opt = obj_simu_anneal(y);
            end
        else
            P = rand(1,1);
            if P < exp(-delta_E/T)
                x = y;
            end
        end
%         ph_iter.XData = x;
%         ph_iter.YData = f(x);
%         te_1.Position = [x+0.3 f(x) 0];
%         legend('show');
%         refresh
%         drawnow;
 
        iter = iter + 1;
    end
    T = T * 0.95;
end