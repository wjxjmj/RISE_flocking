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

% This project consists of the first simulation in the paper.

% If you have questions about this code, please feel free to contact the 
% author of the code or the author of the correspondence.
% E-mail of the author of this code: wjxjmj@126.com

% All parameters are wrapped into the "para" structure.
para=[];
para.kx=1;
para.kv=5;
para.kl=para.kx+para.kv;
para.dim=2;
para.N=50;
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

% Define the structure "state" so that the three systems have the same
% initial position and velocity, thus facilitating comparison.
state=[];
state.x = unifrnd(-3,3,[para.dim,para.N]);
state.v = unifrnd(-1,1,[para.dim,para.N]);
state.s = zeros(size(state.v));

% Configuration for solvers.
solver = 'ode113';
odeopt = odeset('reltol',1e-3,'abstol',1e-6);

% Simulate the proposed algorithm.
sim1 = vectorField(state,@(t,s)rise_flocking(t,s,para));
sim1.set_size_check(true);
sim1.set_options(odeopt);
sim1.solve(solver,para.tspan,state);
sim1.addSignals(@(t,s)sig(t,s,para));
[t1,data1]=sim1.result();


figure(1)
clf;hold on
plot(data1.x(1,1:para.dim:end),data1.x(1,2:para.dim:end),'kx');
plot(data1.x(:,1:para.dim:end),data1.x(:,2:para.dim:end),'c-');
plot(data1.x(end,1:para.dim:end),data1.x(end,2:para.dim:end),'bo', 'MarkerFaceColor','b','MarkerSize',5);
plot(data1.xd(:,1),data1.xd(:,2),'m--');
plot(data1.xd(end,1),data1.xd(end,2),'rp');
xe = reshape(data1.x(end,:),[para.dim,para.N]);
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
export_fig compare_rise_flocking.png -q101 -m4

function signals = sig(t,state,para)
[xd,vd,ad,dad] = leader(t,para);
signals.xd = xd;
signals.vd = vd;
end

function dot_state = rise_flocking(t,state,para)
x = state.x;
v = state.v;
s = state.s;
[xd,vd,ad,dad] = leader(t,para);
[d,dd] = disturbance(t,para);
u = zeros(size(v));
us = zeros(size(s));
E1 = zeros(size(x));
dE1 = zeros(size(x));
E2 = zeros(size(x));

for i=1:para.N
    xi = x(:,i)-xd;
    vi = v(:,i)-vd;
    si = s(:,i);
    E1(:,i) = para.kx * xi; % Eq. (20)
    dE1(:,i)= para.kx * vi;
    for j=1:para.N
        if i==j
            continue
        end
        xj = x(:,j)-xd;
        vj = v(:,j)-vd;
        E1(:,i) = E1(:,i) + phi(xi-xj,50,para.r,para.l); % Eq. (20)
        dE1(:,i)=dE1(:,i) + dPhi(xi-xj,vi-vj,50,para.r,para.l);
    end
    E2(:,i) = para.kl*vi + 1*para.alpha1(i)*E1(:,i); % Eq. (21)
    us(:,i) = -(1+para.ks(i))*para.alpha2(i)*E2(:,i)-para.beta(i)*sign(E2(:,i)); % Eq. (26)
    u(:,i) = -(1+para.ks(i))*E2(:,i)-para.alpha1(i)/para.kl*dE1(:,i)+si; % Eq. (25)
end

dot_state.x = v;
dot_state.v = u - d;
dot_state.s = us;

end


% trajectories of the virtual leader
function [xd,vd,ad,dad]=leader(t,para)
rd = para.rd;
hz = para.hz;
xd=[rd*cos(t*hz);rd*sin(t*hz)];
vd=[-rd*hz*sin(t*hz);rd*hz*cos(t*hz)];
ad=[-rd*hz^2*cos(t*hz);-rd*hz^2*sin(t*hz)];
dad=[rd*hz^3*sin(t*hz);-rd*hz^3*cos(t*hz)];
end

% Eq. (4)
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

% Eq. (6)
% Remark: For ease of implementation, we have modified the potential field 
% functions for ease of derivation, which has very little effect on the 
% conclusions of the paper.
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

% Eq. (3)
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

% bounded disturbances
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


