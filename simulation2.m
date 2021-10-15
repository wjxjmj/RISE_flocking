% Matlab source code for paper "RISE-based Distributed Flocking Control of 
% Double-Integrator Multi-Agent Systems with A Varying Virtual Leader" with
% DOI: XXX

% This project relies on the open source simulation toolkit:
% openVectorField. Before running this code, please ensure that you have
% configured the openVectorField toolkit correctly. This toolkit can be 
% found at https://github.com/wjxjmj/openVectorField.

% This project relies on the open source simulation toolkit:
% export_fig. Before running this code, please ensure that you have
% configured the export_fig toolkit correctly. This toolkit can be 
% found at https://github.com/altmany/export_fig.

% This project consists of three sets of simulations, of which:
% sim1 corresponds to the algorithm proposed in this paper;
% sim2 corresponds to the Ghapani's algorithm;
% sim3 corresponds to the Olfati-Saber's algorithm;

% If you have questions about this code, please feel free to contact the 
% author of the code or the author of the correspondence.
% E-mail of the author of this code: wjxjmj@126.com

% All parameters are wrapped into the "para" structure.
para.kx=1;
para.kv=5;
para.kl=para.kx+para.kv;
para.dim=2;
para.N=8;
para.rd=5;
para.tspan=linspace(0,20,1000);
para.hz=1/20*1.8*pi;

para.a=unifrnd(-1,1,[para.N,1]);
para.b=unifrnd(-1,1,[para.N,1]);

para.r=1.2;
para.l=1;

para.alpha1=ones(para.N,1);
para.alpha2=ones(para.N,1);
para.beta=ones(para.N,1);
para.ks=5*ones(para.N,1);

para.gamma1=ones(para.N,para.N+1);
para.gamma2=ones(1,para.N);

para.alpha1=ones(para.N,1);
para.alpha2=ones(para.N,1);

% Define the structure "state" so that the three systems have the same 
% initial position and velocity, thus facilitating comparison.
state=[];
state.x = unifrnd(-3,3,[para.dim,para.N]);
state.v = unifrnd(-1,1,[para.dim,para.N]);

% Define structure "state1" for the proposed algorithm.
state1=state;
state1.s = zeros(size(state1.v));
state1.P = 10000;
state1.Q = 10000;

% Define structure "state2" for the Ghapani's algorithm.
state2=state;
state2.nu = zeros(size(state2.v));
state2.alpha = ones([para.N,para.N+1]);
state2.beta = ones([1,para.N]);

% Define structure "state3" for the Olfati-Saber's algorithm.
state3 = state;

% Configuration for solvers.
solver = 'ode45';
odeopt = odeset('reltol',1e-3,'abstol',1e-6);
steps = 0.0001;

% Simulate the proposed algorithm.
sim1 = vectorField(state1,@(t,s)rise_flocking(t,s,para));
sim1.set_size_check(true);
% Uncomment this line if "solver" is defined as "Euler".
% sim1.set_dt(steps);
sim1.set_options(odeopt);
sim1.solve(solver,para.tspan,state1);
sim1.addSignals(@(t,s)enegy1(t,s,para));
[t1,data1]=sim1.result();

% Simulate the Ghapani's algorithm.
sim2 = vectorField(state2,@(t,s)fully_flocking(t,s,para));
sim2.set_size_check(true);
% Uncomment this line if "solver" is defined as "Euler".
% sim2.set_dt(steps);
sim2.set_options(odeopt);
sim2.solve(solver,para.tspan,state2);
sim2.addSignals(@(t,s)enegy2(t,s,para));
[t2,data2]=sim2.result();

% Simulate the Olfati-Saber's algorithm.
sim3 = vectorField(state3,@(t,s)olfati_flocking(t,s,para));
sim3.set_size_check(true);
% Uncomment this line if "solver" is defined as "Euler".
% sim3.set_dt(steps);
sim3.set_options(odeopt);
sim3.solve(solver,para.tspan,state3);
sim3.addSignals(@(t,s)enegy3(t,s,para));
[t3,data3]=sim3.result();

% Visulization.
figure(1)
plot(t1,log10(data1.e))
grid on
xlabel('time (s)')
ylabel('$V$','interpreter','latex')
height=151;
set(gcf,'units','points','position',[0,0,height*42/15,height]);
set(gcf, 'Color', 'w');
export_fig compare_rise_lyap.png -q101 -m4

figure(2)
plot(t1,sign(data1.dot_e).*log10(abs(data1.dot_e)))
grid on
xlabel('time (s)')
ylabel('$\dot{V}$','interpreter','latex')
height=151;
set(gcf,'units','points','position',[0,0,height*42/15,height]);
set(gcf, 'Color', 'w');
export_fig compare_rise_dot_lyap.png -q101 -m4

figure(3)
plot(t2,log10(data2.e))
grid on
xlabel('time (s)')
ylabel('$V$','interpreter','latex')
height=151;
set(gcf,'units','points','position',[0,0,height*42/15,height]);
set(gcf, 'Color', 'w');
export_fig compare_fully_lyap.png -q101 -m4

figure(4)
plot(t2,sign(data2.dot_e).*log10(abs(data2.dot_e)))
grid on
xlabel('time (s)')
ylabel('$\dot{V}$','interpreter','latex')
height=151;
set(gcf,'units','points','position',[0,0,height*42/15,height]);
set(gcf, 'Color', 'w');
export_fig compare_fully_dot_lyap.png -q101 -m4

figure(5)
plot(t3,log10(data3.e))
grid on
xlabel('time (s)')
ylabel('$V$','interpreter','latex')
height=151;
set(gcf,'units','points','position',[0,0,height*42/15,height]);
set(gcf, 'Color', 'w');
export_fig compare_olfati_lyap.png -q101 -m4

figure(6)
plot(t3,sign(data3.dot_e).*log10(abs(data3.dot_e)))
grid on
xlabel('time (s)')
ylabel('$\dot{V}$','interpreter','latex')
height=151;
set(gcf,'units','points','position',[0,0,height*42/15,height]);
set(gcf, 'Color', 'w');
export_fig compare_olfati_dot_lyap.png -q101 -m4

figure(7)
clf;hold on
plot(data1.x(1,1:para.dim:end),data1.x(1,2:para.dim:end),'kx');
plot(data1.x(:,1:para.dim:end),data1.x(:,2:para.dim:end),'c-');
plot(data1.x(end,1:para.dim:end),data1.x(end,2:para.dim:end),'bo', 'MarkerFaceColor','b','MarkerSize',5);
plot(data1.xd(:,1),data1.xd(:,2),'m--');
plot(data1.xd(end,1),data1.xd(end,2),'rp');
xe= reshape(data1.x(end,:),[para.dim,para.N]);
for i=1:para.N
        xi = xe(:,i);
        drawCircle(xi,para.l,[0.5,0.6,0.8],'--');
        for j=1:para.N
            if i==j
                continue
            else
                
                xj = xe(:,j);
                dis = norm(xi-xj);
                if dis<=para.l
                    line([xi(1) xj(1)],[xi(2) xj(2)]);
                end
            end
        end
end
hold off
axis equal
grid on
title('flocking of agents')
xlabel('x axis')
ylabel('y axis')
height=400;
set(gcf,'units','points','position',[0,0,height*50/42,height]);
set(gcf, 'Color', 'w');
% export_fig compare_rise_flocking.png -q101 -m4

figure(8)
ndis=getInputs2(data1.u,para);
plot(t1,ndis');
grid on
% title('Time evolution of inputs')
xlabel('seconds')
ylabel('norm of inputs')
height=151;
set(gcf,'units','points','position',[0,0,height*42/15,height]);
set(gcf, 'Color', 'w');
export_fig compare_rise_inputs_ode.png -q101 -m4
% export_fig compare_rise_inputs_small.png -q101 -m4

figure(9)
clf;hold on
plot(data2.x(1,1:para.dim:end),data2.x(1,2:para.dim:end),'kx');
plot(data2.x(:,1:para.dim:end),data2.x(:,2:para.dim:end),'c-');
plot(data2.x(end,1:para.dim:end),data2.x(end,2:para.dim:end),'bo', 'MarkerFaceColor','b','MarkerSize',5);
plot(data2.xd(:,1),data2.xd(:,2),'m--');
plot(data2.xd(end,1),data2.xd(end,2),'rp');
xe= reshape(data2.x(end,:),[para.dim,para.N]);
for i=1:para.N
        xi = xe(:,i);
        drawCircle(xi,para.l,[0.5,0.6,0.8],'--');
        for j=1:para.N
            if i==j
                continue
            else
                
                xj = xe(:,j);
                dis = norm(xi-xj);
                if dis<=para.l %abs(dis-1)<=0.005
                    line([xi(1) xj(1)],[xi(2) xj(2)]);
                end
            end
        end
end
hold off
axis equal
grid on
% title('flocking of agents')
xlabel('x axis')
ylabel('y axis')
height=400;
set(gcf,'units','points','position',[0,0,height*50/42,height]);
set(gcf, 'Color', 'w');
export_fig compare_fully_flocking.png -q101 -m4

figure(10)
ndis=getInputs2(data2.u,para);
plot(t2,ndis');
grid on
% title('Time evolution of inputs')
xlabel('seconds')
ylabel('norm of inputs')
height=151;
set(gcf,'units','points','position',[0,0,height*42/15,height]);
set(gcf, 'Color', 'w');
export_fig compare_fully_inputs_ode.png -q101 -m4

figure(11)
clf;hold on
plot(data3.x(1,1:para.dim:end),data3.x(1,2:para.dim:end),'kx');
plot(data3.x(:,1:para.dim:end),data3.x(:,2:para.dim:end),'c-');
plot(data3.x(end,1:para.dim:end),data3.x(end,2:para.dim:end),'bo', 'MarkerFaceColor','b','MarkerSize',5);
plot(data3.xd(:,1),data3.xd(:,2),'m--');
plot(data3.xd(end,1),data3.xd(end,2),'rp');
xe= reshape(data3.x(end,:),[para.dim,para.N]);
for i=1:para.N
        xi = xe(:,i);
        drawCircle(xi,para.l,[0.5,0.6,0.8],'--');
        for j=1:para.N
            if i==j
                continue
            else
                
                xj = xe(:,j);
                dis = norm(xi-xj);
                if dis<=para.l
                    line([xi(1) xj(1)],[xi(2) xj(2)]);
                end
            end
        end
end
hold off
axis equal
grid on
% title('flocking of agents')
xlabel('x axis')
ylabel('y axis')
height=400;
set(gcf,'units','points','position',[0,0,height*50/42,height]);
set(gcf, 'Color', 'w');
export_fig compare_olfati_flocking.png -q101 -m4

figure(12)
ndis=getInputs2(data3.u,para);
plot(t3,ndis');
grid on
% title('Time evolution of inputs')
xlabel('seconds')
ylabel('norm of inputs')
height=151;
set(gcf,'units','points','position',[0,0,height*42/15,height]);
set(gcf, 'Color', 'w');
export_fig compare_olfati_inputs_ode.png -q101 -m4
% export_fig compare_olfati_inputs_small.png -q101 -m4





function signals = enegy1(t,state,para)
[~ , signals] = rise_flocking(t,state,para);
end

function [dot_state , signals] = rise_flocking(t,state,para)
x = state.x;
v = state.v;
s = state.s;
P = state.P;
Q = state.Q;
[xd,vd,ad,dad] = leader(t,para);
[d,dd] = disturbance(t,para);
u = zeros(size(v));
du = zeros(size(v));
us = zeros(size(s));
E1 = zeros(size(x));
dE1 = zeros(size(x));
ddE1 = zeros(size(x));
E2 = zeros(size(x));
dE2 = zeros(size(x));
E3 = zeros(size(x));
dE3 = zeros(size(x));
s1=zeros([1,para.N]);
s2=zeros([1,para.N]);
beta_sign_e2=zeros(size(s));
phi_ = zeros(size(x));
DPhi = zeros(size(x));
Psi=zeros([1,para.N]);

for i=1:para.N
    xi = x(:,i)-xd;
    vi = v(:,i)-vd;
    si = s(:,i);
    E1(:,i) = para.kx * xi;
    dE1(:,i)= para.kx * vi;
    for j=1:para.N
        if i==j
            continue
        end
        xj = x(:,j)-xd;
        vj = v(:,j)-vd;
        E1(:,i) = E1(:,i) + phi(xi-xj,50,para.r,para.l);
        dE1(:,i)=dE1(:,i) + dPhi(xi-xj,vi-vj,50,para.r,para.l);
        phi_(:,i)=phi_(:,i) + phi(xi-xj,50,para.r,para.l);
        DPhi(:,i)=DPhi(:,i) + dPhi(xi-xj,vi-vj,50,para.r,para.l);
        Psi(i) = Psi(i) + psi(norm(xi-xj)^2,50,para.r,para.l);
    end
    E2(:,i) = para.kl*vi + 1*para.alpha1(i)*E1(:,i);
    
    us(:,i) = -(1+para.ks(i))*para.alpha2(i)*E2(:,i)-para.beta(i)*sign(E2(:,i));
    u(:,i) = -(1+para.ks(i))*E2(:,i)-para.alpha1(i)/para.kl*dE1(:,i)+si;
    dE2(:,i) = para.kl*(u(:,i)-ad-d(:,i)) + 1*para.alpha1(i)*dE1(:,i);
    beta_sign_e2(:,i) = para.beta(i)*sign(E2(:,i));
end

for i=1:para.N
    xi = x(:,i)-xd;
    vi = v(:,i)-vd;
    ui = u(:,i);
    ddE1(:,i) = para.kx * ui;
    si = s(:,i);
    for j=1:para.N
        if i==j
            continue
        end
        xj = x(:,j)-xd;
        vj = v(:,j)-vd;
        uj = u(:,j);
        ddE1(:,i)=ddE1(:,i) + ddPhi(xi-xj,vi-vj,ui-uj,50,para.r,para.l);
    end
    
    du(:,i) = -(1+para.ks(i))*dE2(:,i)-para.alpha1(i)/para.kl*ddE1(:,i)+us(:,i);
    
    E3(:,i) = dE2(:,i) + para.alpha2(i)*E2(:,i);
    
    dE3(:,i) = para.kl*(du(:,i)-dad-dd(:,i)) + 1*para.alpha1(i)*ddE1(:,i)...
               + para.alpha2(i)*dE2(:,i);
    
end

u_P=0;
u_Q=0;
for i=1:para.N
    hd = -para.kl * dad - para.kl* dd(:,i);
    u_P = u_P -E3(:,i)'*(hd-para.kl*beta_sign_e2(:,i)); 
    u_Q = u_Q -E1(:,i)'*DPhi(:,i);
end

e=0;
dot_e=0;

k1=1;
k2=1;
k3=1;
k4=1;
k5=1;
k6=1;
k7=1;

for i=1:para.N
    xi = x(:,i)-xd;
    vi = v(:,i)-vd;
    e = e + k1 * 0.5*norm(E1(:,i))^2 ...
          + k2 * 0.5*norm(E2(:,i))^2 ...
          + k3 * 0.5*norm(E3(:,i))^2 ...
          + k4 * para.kv*Psi(i)...
          + k5 * para.kx*0.5*(xi'*xi);
    dot_e = dot_e + k1 *E1(:,i)'*(E2(:,i)-para.alpha1(i)*E1(:,i)+DPhi(:,i)-para.kv*(vi)) ...
                  + k2 *E2(:,i)'*( E3(:,i) - para.alpha2(i)*E2(:,i))...
                  + k3 *E3(:,i)'*(-para.kl*(1+para.ks(i))*E3(:,i))...
                  + k3 *E3(:,i)'*(-para.kl*beta_sign_e2(:,i))...
                  + k3 *E3(:,i)'*(-E1(:,i)-E2(:,i))...
                  + k3 *E3(:,i)'*(para.alpha2(i)*E3(:,i))...
                  + k3 *E3(:,i)'*((1 - para.alpha2(i)^2)*E2(:,i)+E1(:,i))...
                  + k5 *para.kx*(xi'*vi)...
                  - k6 *E3(:,i)'*(-para.kl*(dad+dd(:,i)))...
                  + k3 *E3(:,i)'*(-para.kl*(dad+dd(:,i)))...
                  - k6 *E3(:,i)'*(-para.kl*beta_sign_e2(:,i))...
                  - k7 *E1(:,i)'*DPhi(:,i);
              for j=1:para.N
                  xj = x(:,j)-xd;
                  vj = v(:,j)-vd;
                  if i==j
                      continue
                  else
                      xij = xi-xj;
                      vij = vi-vj;
                      dot_e = dot_e + k4 * para.kv* phi(xij,50,para.r,para.l)'*vij;
                  end
              end
end

e = e +k6*P+k7*Q;


de1 = E1(:)'*dE1(:);
de2 = E2(:)'*dE2(:);
de3 = E3(:)'*dE3(:);

signals.e = e;
signals.dot_e =dot_e;
signals.e1 = 0.5*norm(E1(:))^2;
signals.e2 = 0.5*norm(E2(:))^2;
signals.e3 = 0.5*norm(E3(:))^2;
signals.P = P;
signals.Q = Q;
signals.Psi = Psi;
signals.xd = xd;
signals.u = u;
signals.d=d;
signals.E1 = E1;
signals.E2 = E2;
signals.E3 = E3;
signals.s1 = s1;
signals.dE1 = dE1;
signals.dE2 = dE2;
signals.dE3 = dE3;
signals.s2 = s2;
signals.ddE1 = ddE1;

dot_state.x = v;
dot_state.v = u - d;
dot_state.s = us;




dot_state.P = u_P;
dot_state.Q = u_Q;

end

function signals = enegy2(t,state,para)
[~ , signals] = fully_flocking(t,state,para);
end

function [dot_state, signals] = fully_flocking(t,state,para)
x = state.x;
v = state.v;
nu = state.nu;
alpha = state.alpha;
beta = state.beta;

u_alpha = zeros(size(alpha));
u_beta = zeros(size(beta));
u_nu = zeros(size(nu));

[xd,vd,ad,dad] = leader(t,para);
tao_l = norm(ad);
[d,dd] = disturbance(t,para);
u = zeros(size(v));
V1 = 0;
dV1 = 0;
V2 = 0;
dV2 = 0;
V3 = 0;
dV3 = 0;

action = zeros(size(x));
consensus = zeros(size(x));
sliding = zeros(size(x));

bar_alpha=50;
bar_beta=tao_l+1;

for i=1:para.N
    xi = x(:,i);
    vi = v(:,i);
    nui = nu(:,i);
    bi=1;
    V2 = 0.5 * (xi'*xi);
    dV2= xi'*vi;
    for j=1:para.N
        if i==j
            continue
        end
        xj = x(:,j);
        vj = v(:,j);
        dis2 = (xi-xj)'*(xi-xj);
        aij = rho(dis2/para.r^2);
        V2 = V2 + psi(norm(xi-xj)^2,50,para.r,para.l);
        dV2=dV2 + phi(xi-xj,50,para.r,para.l)'*(vi-vj);
        action(:,i)=action(:,i) + phi(xi-xj,50,para.r,para.l);
        consensus(:,i)=consensus(:,i)+alpha(i,j)*aij*sign(vi-vj);
        
        u_alpha(i,j)=para.gamma1(i)*aij*norm(vi-vd,1);
        V3=V3+1/4/para.gamma1(i)*(alpha(i,j)-bar_alpha)^2;
        dV3=dV3+1/2/para.gamma1(i)*(alpha(i,j)-bar_alpha)*u_alpha(i,j);

    end
    si = vi - nui;
    
    sliding(:,i)=beta(i)*sign(si);
    u_beta(i)=para.gamma2(i)*norm(si,1);
    u_alpha(i,para.N+1)=para.gamma1(i)*bi*norm(vi-vd,1);
    
    u_nu(:,i) = -action(:,i) - consensus(:,i)...
                -bi*(xi-xd)-bi*alpha(i,para.N+1)*sign(vi-vd);
    u(:,i)    = -action(:,i) - consensus(:,i)...
                -bi*(xi-xd)-bi*alpha(i,para.N+1)*sign(vi-vd)...
                -sliding(:,i);
            
    V1=V1 + 0.5* (si'*si);
    dV1=dV1+si'*(u(:,i)-d(:,i)-u_nu(:,i));
    V3=V3+1/2/para.gamma1(i)*(alpha(i,para.N+1)-bar_alpha)^2;
    V3=V3+1/2/para.gamma2(i)*(beta(i)-bar_beta)^2;
    
    dV3=dV3+1/para.gamma1(i)*(alpha(i,para.N+1)-bar_alpha)'*u_alpha(i,para.N+1);
    dV3=dV3+1/para.gamma2(i)*(beta(i)-bar_beta)'*u_beta(i);
    
end

signals.e = V1+V2+V3;
signals.dot_e =dV1+dV2+dV3;

signals.xd = xd;
signals.u = u;

dot_state.x = v;
dot_state.v = u - d;
dot_state.nu = u_nu;
dot_state.alpha = u_alpha;
dot_state.beta = u_beta;


end

function signals = enegy3(t,state,para)
[~ , signals] = olfati_flocking(t,state,para);
end

function [dot_state, signals] = olfati_flocking(t,state,para)
x = state.x;
v = state.v;


[xd,vd,ad,dad] = leader(t,para);
[d,dd] = disturbance(t,para);
u = zeros(size(v));

action = zeros(size(x));
consensus = zeros(size(x));
Psi=zeros([1,para.N]);
V1=zeros([1,para.N]);
V2=zeros([1,para.N]);
dV1=zeros([1,para.N]);
dV2=zeros([1,para.N]);
Lyap=0;
dot_Lyap=0;
for i=1:para.N
    xi = x(:,i);
    vi = v(:,i);
    V1(i) = 0.5*norm(xi-xd)^2;
    V2(i) = 0.5*norm(vi-vd)^2;
    
    for j=1:para.N
        if i==j
            continue
        end
        xj = x(:,j);
        vj = v(:,j);
        dis2 = (xi-xj)'*(xi-xj);
        aij = rho(dis2/para.r^2);
        action(:,i)=action(:,i) + phi(xi-xj,50,para.r,para.l);
        consensus(:,i)=consensus(:,i)+aij*(vi-vj);
        
        Psi(i) = Psi(i) + psi(norm(xi-xj)^2,50,para.r,para.l);

    end
    
    u(:,i) = -action(:,i) - consensus(:,i) - (xi-xd+vi-vd) +ad;
    ui = u(:,i);
    
    dV1(:,i)= (xi-xd)'*(vi-vd);
    dV2(:,i)= (vi-vd)'*(ui-ad);
    
    k1=1;
    k2=1;
    k3=1;
    
    Lyap = Lyap + k1* V1(i) + k2* V2(i) + k3* Psi(i);
    dot_Lyap = dot_Lyap + k1* dV1(:,i) + k2* dV2(:,i);
    for j=1:para.N
        xj = x(:,j);
        vj = v(:,j);
        if i==j
            continue
        else
            xij = xi-xj;
            vij = vi-vj;
            dot_Lyap = dot_Lyap + k3* phi(xij,50,para.r,para.l)'*vij;
        end
    end
end

signals.e = Lyap;
signals.dot_e =dot_Lyap;

signals.V1 = V1;
signals.V2 = V2;
signals.Psi = Psi;
signals.xd = xd;
signals.d=d;
signals.u = u;


dot_state.x = v;
dot_state.v = u - d;



end

function [xd,vd,ad,dad]=leader(t,para)
rd = para.rd;
hz = para.hz;
xd=[rd*cos(t*hz);rd*sin(t*hz)];
vd=[-rd*hz*sin(t*hz);rd*hz*cos(t*hz)];
ad=[-rd*hz^2*cos(t*hz);-rd*hz^2*sin(t*hz)];
dad=[rd*hz^3*sin(t*hz);-rd*hz^3*cos(t*hz)];
end

function y = psi(x,k,r,l)
if x<r^2
    y =0.5* (1/24).*k.*pi.^(-3).*r.^(-2).*(4.*pi.^3.*x.*((3.*r.^2+(-2).*x).*x+ ...
  3.*l.^2.*((-2).*r.^2+x))+(-6).*pi.*r.^4.*((-1).*l.^2+x).*cos(2.* ...
  pi.*r.^(-2).*x)+3.*r.^6.*sin(2.*pi.*r.^(-2).*x));
else
    y=0.5*(1/12).*k.*pi.^(-2).*r.^2.*(l.^2.*(3+(-6).*pi.^2)+((-3)+2.*pi.^2) ...
  .*r.^2);
end
y=y+(-1/24).*k.*(4.*l.^4.*((-3)+l.^2.*r.^(-2))+3.*pi.^(-3).*r.^4.*sin( ...
  2.*l.^2.*pi.*r.^(-2)));
end



function y = phi(x,k,r,l)
dis2 = x'*x;
y = k * rho(dis2/r^2)*(dis2-l^2)*x;
end

function y = dPhi(x,v,k,r,l)
dis2 = x'*x;
y =  k * dRho(dis2/r^2)*2*x'*v/r^2*(dis2-l^2)*x+...
     k * rho(dis2/r^2)*(2*x'*v)*x+...
     k * rho(dis2/r^2)*(dis2-l^2)*v;
end

function y = ddPhi(x,v,u,k,r,l)
dis2 = x'*x;
y =  k * ddRho(dis2/r^2)*(2*x'*v/r^2)^2*(dis2-l^2)*x+...
     k * dRho(dis2/r^2)*2*(v'*v)/r^2*(dis2-l^2)*x+...
     k * dRho(dis2/r^2)*2*x'*u/r^2*(dis2-l^2)*x+...
     k * dRho(dis2/r^2)*2*x'*v/r^2*(2*x'*v-l^2)*x+...
     k * dRho(dis2/r^2)*2*x'*v/r^2*(dis2-l^2)*v+...
     k * dRho(dis2/r^2)*(2*x'*v/r^2)*(2*x'*v)*x+...
     k * rho(dis2/r^2)*2*(v'*v)*x+...
     k * rho(dis2/r^2)*(2*x'*u)*x+...
     k * rho(dis2/r^2)*(2*x'*v)*v+...
     k * dRho(dis2/r^2)*(2*x'*v/r^2)*(dis2-l^2)*v+...
     k * rho(dis2/r^2)*(2*x'*v-l^2)*v+...
     k * rho(dis2/r^2)*(dis2-l^2)*u;
end

function y = rho(x)
if x<1
    y=1-x+sin(2*pi*x)/2/pi;
else
    y=0;
end
end

function y = dRho(x)
if x<1
    y=-1 + cos(2*pi*x);
else
    y=0;
end
end

function y = ddRho(x)
if x<1
    y= -2*pi * sin(2*pi*x);
else
    y=0;
end
end

function [d,dd] = disturbance(t,para)
a = para.a;
b = para.b;
d = zeros(para.dim,para.N);
dd = d;

for i=1:para.N
    d(:,i) = 5*[sin(t+a(i)*pi);cos(t+b(i)*pi)];
    dd(:,i)= 5*[cos(t+a(i)*pi);-sin(t+b(i)*pi)];
end
end

function y=getInputs2(uHis,para)
[loop,~]=size(uHis);

dis=zeros(para.N,loop);
iter=1;
for t=1:loop
    u=reshape(uHis(t,:),[para.dim,para.N]);
    for i=1:para.N
        dis(i,iter)=norm(u(:,i));
    end
    iter=iter+1;
end
y=dis;
end


