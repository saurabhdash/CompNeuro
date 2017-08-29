function Project2_MLE
close all;
clear;
clc;
%% MLE Constants
global gca;
gca = 4.4;
global gk
gk = 8;
global gl
gl = 2;
global Vca
Vca = 120;
global Vk
Vk = -84;
global Vl
Vl = -60;
global phi
phi = 0.02;
global V1
V1 = -1.2;
global V2
V2 = 18;
global V3
V3 = 2;
global V4
V4 = 30;
global V5
V5 = 2;
global V6
V6 = 30;
global C
C = 20;
global Iext;
Iext = 0;
%%
syms V w m_inf w_inf Tw;
m_inf = 0.5*(1+tanh((V-V1)/V2));
w_inf = 0.5*(1+tanh((V-V3)/V4));
Tw = 1/cosh((V-V3)/(2*V4));
%% Part 2
eqns = [ phi*(w_inf-w)/Tw == 0, gca*m_inf*(V-Vca)+ gk*w*(V-Vk) + gl*(V-Vl) == 0];
vars = [V w];
S = vpasolve(eqns, vars);
[S.V S.w];
disp(['For Iext = ' num2str(Iext) ', the Equilibrium point is V = ' num2str(double(S.V)) 'mV & w = ' num2str(double(S.w))])
%%
w1 = (-gca*m_inf*(V-Vca)-gl*(V-Vl))/(gk*(V-Vk));
V = -83:0.1:120;
w_nc = double(subs(w_inf));
V_nc = double(subs(w1));
figure(4)
plot(V,100*w_nc); % w Null cline
hold on;
plot(V,100*V_nc); % V Null cline
xlim([-84 120])
ylim([0 100])
title('Null Clines and Quiver Plot')
xlabel('V(mV)') % x-axis label
ylabel('w') % y-axis label
legend('w Null Cline','V Null Cline')

figure(1)
plot(V,w_nc); % w Null cline
hold on;
plot(V,V_nc); % V Null cline
xlim([-84 120])
ylim([0 1])
title('Null Clines')
xlabel('V(mV)') % x-axis label
ylabel('w') % y-axis label
legend('w Null Cline','V Null Cline')
%% Part 3
syms V;
J = jacobian([(1/C)*(Iext-(gca*m_inf*(V-Vca)+ gk*w*(V-Vk) + gl*(V-Vl))), phi*(w_inf-w)/Tw], [V, w]);
V = S.V;
w = S.w;
Jeq = double(subs(J));
Eign_eq = eig(Jeq);
disp(['For V = ' num2str(double(V)) ' & w = ' num2str(double(w)) ', the Eigenvalues are = ' num2str(Eign_eq(1)) ' & '  num2str(Eign_eq(2))])
%%
[V,w] = meshgrid(-84:10:120, 0:0.1:1);
m_inf = 0.5*(1+tanh((V-V1)./V2));
w_inf = 0.5*(1+tanh((V-V3)./V4));
Tw = 1./cosh((V-V3)./(2*V4));
dw = phi.*(w_inf-w)./Tw;
dv = (1/C)*(Iext-gca.*m_inf.*(V-Vca) - gk.*w.*(V-Vk) - gl.*(V-Vl));
figure(4);
hold on;
quiver(V,100*w,dv,100*dw)
%% Part 5 and 6
syms V w m_inf w_inf Tw;
m_inf = 0.5*(1+tanh((V-V1)/V2));
w_inf = 0.5*(1+tanh((V-V3)/V4));
Tw = 1/cosh((V-V3)/(2*V4));
[t1,S1]=ode15s(@MLE, [0 300], [0  , 0.014]);
[t2,S2]=ode15s(@MLE, [0 300], [-10, 0.014]);
[t3,S3]=ode15s(@MLE, [0 300], [-25, 0.014]);
figure(1)
plot(S1(:,1), S1(:,2),'b',S2(:,1), S2(:,2),'b',S3(:,1), S3(:,2),'b');
title('Phase Plots')
xlabel('V(mV)') % x-axis label
ylabel('w') % y-axis label
%legend('w Null Cline','V Null Cline')
figure;
plot(t1,S1(:,1),'b');
hold on;
plot(t2,S2(:,1),'b');
plot(t3,S3(:,1),'b');
title('Transient Response')
xlabel('t(ms)') % x-axis label
ylabel('V(mV)') % y-axis label
%legend(['Initial Voltage = ' num2str(0)],['Initial Voltage = ' num2str(-10)],['Initial Voltage = ' num2str(-25)])

phi = 0.04;
[t1,S1]=ode15s(@MLE, [0 300], [0  , 0.014]);
[t2,S2]=ode15s(@MLE, [0 300], [-10, 0.014]);
[t3,S3]=ode15s(@MLE, [0 300], [-25, 0.014]);
figure(1);
plot(S1(:,1), S1(:,2),'r',S2(:,1), S2(:,2),'r',S3(:,1), S3(:,2),'r');
figure(2);
plot(t1,S1(:,1),'r');
hold on;
plot(t2,S2(:,1),'r');
plot(t3,S3(:,1),'r');

phi = 0.01;
[t1,S1]=ode15s(@MLE, [0 300], [0  , 0.014]);
[t2,S2]=ode15s(@MLE, [0 300], [-10, 0.014]);
[t3,S3]=ode15s(@MLE, [0 300], [-25, 0.014]);
figure(1);
plot(S1(:,1), S1(:,2),'g',S2(:,1), S2(:,2),'g',S3(:,1), S3(:,2),'g');
figure(2);
plot(t1,S1(:,1),'g');
hold on;
plot(t2,S2(:,1),'g');
plot(t3,S3(:,1),'g');
%% Part 6
figure;
phi = 0.02;
max_v = 0;
for v = -17:0.1:-10
    [t1,S1]=ode15s(@MLE, [0 300], [v  , 0.014]);
    %plot(S1(:,1), S1(:,2));
    %hold on;
    max_v = [max_v max(S1(:,1))];
end
plot(-17:0.1:-10,max_v(2:end));
title('Threshold Plot')
xlabel('V(mV)') % x-axis label
ylabel('V_max (mV)') % y-axis label
%%
Iext = 86;
eqns = [ phi*(w_inf-w)/Tw == 0, Iext - gca*m_inf*(V-Vca) - gk*w*(V-Vk) - gl*(V-Vl) == 0];
vars = [V w];
S = vpasolve(eqns, vars);
[S.V S.w];
disp(['For Iext = ' num2str(Iext) ', the Equilibrium point is V = ' num2str(double(S.V)) 'mV & w = ' num2str(double(S.w))])
w1 = (Iext-gca*m_inf*(V-Vca)-gl*(V-Vl))/(gk*(V-Vk));
V = -83:0.1:120;
w_nc = double(subs(w_inf));
V_nc = double(subs(w1));
figure;
plot(V,w_nc); % w Null cline
hold on;
plot(V,V_nc); % V Null cline
xlim([-84 120])
ylim([0 1])
figure(6);
plot(V,w_nc); % w Null cline
hold on;
plot(V,V_nc); % V Null cline
xlim([-84 120])
ylim([0 1])
%%
syms V;
J = jacobian([(1/C)*(Iext-(gca*m_inf*(V-Vca)+ gk*w*(V-Vk) + gl*(V-Vl))), phi*(w_inf-w)/Tw], [V, w]);
V = S.V;
w = S.w;
Jeq = double(subs(J));
Eign_eq = eig(Jeq);
disp(['For V = ' num2str(double(V)) ' & w = ' num2str(double(w)) ', the Eigenvalues are = ' num2str(Eign_eq(1)) ' & '  num2str(Eign_eq(2))])
%%
[t1,S1]=ode15s(@MLE, [0 300], [ -60.855, 0.014]);
[t2,S2]=ode15s(@MLE, [0 300], [ -27.952, 0.119]);
[t3,S3]=ode15s(@MLE, [0 300], [-27.9, 0.17]);
figure(5);
hold on;
plot(S1(:,1), S1(:,2),'r',S2(:,1), S2(:,2),'b',S3(:,1), S3(:,2),'g');
figure;
plot(t1,S1(:,1),'r');
hold on;
plot(t2,S2(:,1),'b');
plot(t3,S3(:,1),'g');
figure(5);
hold on;
[t3,S3]=ode15s(@MLE, [2000 0], [-27.952, 0.119]);
K_ms = convhull(S3(:,1),S3(:,2));
plot(S3(K_ms,1), S3(K_ms,2));
figure(6);
hold on;
plot(S3(:,1), S3(:,2));
%hold on;
%[t3,S3]=ode15s(@MLE, [0 1000], [-20, 0.1376]);
%plot(S3(:,1), S3(:,2));
%[t3,S3]=ode15s(@MLE, [0 1000], [-19.8, 0.1376]);
%plot(S3(:,1), S3(:,2));

max_v = 0;
%figure;
for v = -22:0.1:-17
    [t1,S1]=ode15s(@MLE, [0 1000], [v  , 0.1376]);
    %plot(S1(:,1), S1(:,2));
    %xlim([-84 120])
    %ylim([0 1])
    %hold on;
    max_v = [max_v max(S1(:,1))];
end
figure;
plot(-22:0.1:-17,max_v(2:end));
%% Part 9
syms V w m_inf w_inf Tw;
m_inf = 0.5*(1+tanh((V-V1)/V2));
w_inf = 0.5*(1+tanh((V-V3)/V4));
Tw = 1/cosh((V-V3)/(2*V4));

Iext = 80;
eqns = [ phi*(w_inf-w)/Tw == 0, Iext - gca*m_inf*(V-Vca) - gk*w*(V-Vk) - gl*(V-Vl) == 0];
vars = [V w];
S = vpasolve(eqns, vars);
[S.V S.w];
disp(['For Iext = ' num2str(Iext) ', the Equilibrium point is V = ' num2str(double(S.V)) 'mV & w = ' num2str(double(S.w))])
w1 = (Iext-gca*m_inf*(V-Vca)-gl*(V-Vl))/(gk*(V-Vk));
V = -83:0.1:120;
w_nc = double(subs(w_inf));
V_nc = double(subs(w1));
figure;
plot(V,w_nc); % w Null cline
hold on;
plot(V,V_nc); % V Null cline
xlim([-84 120])
ylim([0 1])
syms V;
J = jacobian([(1/C)*(Iext-(gca*m_inf*(V-Vca)+ gk*w*(V-Vk) + gl*(V-Vl))), phi*(w_inf-w)/Tw], [V, w]);
V = S.V;
w = S.w;
Jeq = double(subs(J));
Eign_eq = eig(Jeq);
disp(['For V = ' num2str(double(V)) ' & w = ' num2str(double(w)) ', the Eigenvalues are = ' num2str(Eign_eq(1)) ' & '  num2str(Eign_eq(2))])

syms V w
Iext = 86;
eqns = [ phi*(w_inf-w)/Tw == 0, Iext - gca*m_inf*(V-Vca) - gk*w*(V-Vk) - gl*(V-Vl) == 0];
vars = [V w];
S = vpasolve(eqns, vars);
[S.V S.w];
disp(['For Iext = ' num2str(Iext) ', the Equilibrium point is V = ' num2str(double(S.V)) 'mV & w = ' num2str(double(S.w))])
w1 = (Iext-gca*m_inf*(V-Vca)-gl*(V-Vl))/(gk*(V-Vk));
V = -83:0.1:120;
w_nc = double(subs(w_inf));
V_nc = double(subs(w1));
figure;
plot(V,w_nc); % w Null cline
hold on;
plot(V,V_nc); % V Null cline
xlim([-84 120])
ylim([0 1])
syms V;
J = jacobian([(1/C)*(Iext-(gca*m_inf*(V-Vca)+ gk*w*(V-Vk) + gl*(V-Vl))), phi*(w_inf-w)/Tw], [V, w]);
V = S.V;
w = S.w;
Jeq = double(subs(J));
Eign_eq = eig(Jeq);
disp(['For V = ' num2str(double(V)) ' & w = ' num2str(double(w)) ', the Eigenvalues are = ' num2str(Eign_eq(1)) ' & '  num2str(Eign_eq(2))])

syms V w
Iext = 90;
eqns = [ phi*(w_inf-w)/Tw == 0, Iext - gca*m_inf*(V-Vca) - gk*w*(V-Vk) - gl*(V-Vl) == 0];
vars = [V w];
S = vpasolve(eqns, vars);
[S.V S.w];
disp(['For Iext = ' num2str(Iext) ', the Equilibrium point is V = ' num2str(double(S.V)) 'mV & w = ' num2str(double(S.w))])
w1 = (Iext-gca*m_inf*(V-Vca)-gl*(V-Vl))/(gk*(V-Vk));
V = -83:0.1:120;
w_nc = double(subs(w_inf));
V_nc = double(subs(w1));
figure;
plot(V,w_nc); % w Null cline
hold on;
plot(V,V_nc); % V Null cline
xlim([-84 120])
ylim([0 1])
syms V;
J = jacobian([(1/C)*(Iext-(gca*m_inf*(V-Vca)+ gk*w*(V-Vk) + gl*(V-Vl))), phi*(w_inf-w)/Tw], [V, w]);
V = S.V;
w = S.w;
Jeq = double(subs(J));
Eign_eq = eig(Jeq);
disp(['For V = ' num2str(double(V)) ' & w = ' num2str(double(w)) ', the Eigenvalues are = ' num2str(Eign_eq(1)) ' & '  num2str(Eign_eq(2))])

%% Part 10
gca = 4;
gk = 8;
gl = 2;
Vca = 120;
Vk = -84;
Vl = -60;
phi = 0.0667;
V1 = -1.2;
V2 = 18;
V3 = 12;
V4 = 17.4;
C = 20;
Iext = 30;

syms V w m_inf w_inf Tw;
m_inf = 0.5*(1+tanh((V-V1)/V2));
w_inf = 0.5*(1+tanh((V-V3)/V4));
Tw = 1/cosh((V-V3)/(2*V4));

eqns = [ phi*(w_inf-w)/Tw == 0, Iext - gca*m_inf*(V-Vca) - gk*w*(V-Vk) - gl*(V-Vl) == 0];
vars = [V w];
A = [];
for n = 1:20 
  S = vpasolve(eqns, vars,'random',true); 
  A = [A; double(S.V) double(S.w)];
end
B = unique(A,'rows');
disp(['For Iext = ' num2str(Iext) ', the Equilibrium point is V = ' num2str(B(1,1)) 'mV & w = ' num2str(B(1,2))])
syms V w;
J = jacobian([(1/C)*(Iext-(gca*m_inf*(V-Vca)+ gk*w*(V-Vk) + gl*(V-Vl))), phi*(w_inf-w)/Tw], [V, w]);
V = B(1,1);
w = B(1,2);
Jeq = double(subs(J));
Eign_eq = eig(Jeq);
disp(['For V = ' num2str(double(B(1,1))) ' & w = ' num2str(double(B(1,2))) ', the Eigenvalues are = ' num2str(Eign_eq(1)) ' & '  num2str(Eign_eq(2))])
disp(['For Iext = ' num2str(Iext) ', the Equilibrium point is V = ' num2str(B(2,1)) 'mV & w = ' num2str(B(2,2))])
syms V w;
J = jacobian([(1/C)*(Iext-(gca*m_inf*(V-Vca)+ gk*w*(V-Vk) + gl*(V-Vl))), phi*(w_inf-w)/Tw], [V, w]);
V = B(2,1);
w = B(2,2);
Jeq = double(subs(J));
Eign_eq = eig(Jeq);
disp(['For V = ' num2str(double(B(2,1))) ' & w = ' num2str(double(B(2,2))) ', the Eigenvalues are = ' num2str(Eign_eq(1)) ' & '  num2str(Eign_eq(2))])
disp(['For Iext = ' num2str(Iext) ', the Equilibrium point is V = ' num2str(B(3,1)) 'mV & w = ' num2str(B(3,2))])
syms V w;
J = jacobian([(1/C)*(Iext-(gca*m_inf*(V-Vca)+ gk*w*(V-Vk) + gl*(V-Vl))), phi*(w_inf-w)/Tw], [V, w]);
V = B(3,1);
w = B(3,2);
Jeq = double(subs(J));
Eign_eq = eig(Jeq);
disp(['For V = ' num2str(double(B(3,1))) ' & w = ' num2str(double(B(3,2))) ', the Eigenvalues are = ' num2str(Eign_eq(1)) ' & '  num2str(Eign_eq(2))])
syms V w;
w1 = (Iext-gca*m_inf*(V-Vca)-gl*(V-Vl))/(gk*(V-Vk));
V = -83:0.1:120;
w_nc = double(subs(w_inf));
V_nc = double(subs(w1));
figure;
plot(V,w_nc); % w Null cline
hold on;
plot(V,V_nc); % V Null cline
xlim([-84 120])
ylim([0 1])
syms V w;
[t3,S3]=ode15s(@MLE, [2000 0], [-40.8452, 0.0020475]);
plot(S3(:,1), S3(:,2),'b');
[t3,S3]=ode15s(@MLE, [2000 0], [-40.8452, 0.003]);
plot(S3(:,1), S3(:,2),'b');
[t3,S3]=ode15s(@MLE, [2000 0], [-42.8452, 0.003]);
plot(S3(:,1), S3(:,2),'b');
[t3,S3]=ode15s(@MLE, [2000 0], [-42.8452, 0.001]);
plot(S3(:,1), S3(:,2),'b');

% [t3,S3]=ode15s(@MLE, [2000 0], [-19.6, 0.025]);
% plot(S3(:,1), S3(:,2),'r');
% [t3,S3]=ode15s(@MLE, [2000 0], [-19.5, 0.025]);
% plot(S3(:,1), S3(:,2),'r');
% [t3,S3]=ode15s(@MLE, [2000 0], [-19.5, 0.026]);
% plot(S3(:,1), S3(:,2),'r');
% [t3,S3]=ode15s(@MLE, [0 2000], [-19.62, 0.0255]);
% plot(S3(:,1), S3(:,2),'r');


  end

function dS = MLE(t,S)
global gca;
global gk;
global gl;
global Vca;
global Vk;
global Vl;
global phi;
global V1;
global V2;
global V3;
global V4;
global V5;
global V6;
global C;
global Iext;

V=S(1);
w=S(2);
m_inf = 0.5*(1+tanh((V-V1)/V2));
w_inf = 0.5*(1+tanh((V-V3)/V4));
Tw = 1/cosh((V-V3)/(2*V4));
 
ddt_V = (Iext-(gca*m_inf*(V-Vca)+ gk*w*(V-Vk) + gl*(V-Vl)))/C;
ddt_w = phi*(w_inf-w)/Tw;
dS=[ddt_V; ddt_w];
end