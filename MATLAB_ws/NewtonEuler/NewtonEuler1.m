clc;
clear all;
close all;
clf;
ts = 0.001;
time_span = 0:ts:10;

qc = [sin(time_span)' cos(2*time_span)' sin(3*time_span)'];
qcdot = gradient(qc', ts)';
qcddot = gradient(qcdot', ts)';
% complete inverse dynamics
tau = NewtonEuler(qc, qcdot, qcddot);
% tau = complete dynamics, tau1= gravity, tau2= coriolis, tau_inertia=
% inertia term
% gravity terms
qcdot1=zeros(size(qc));
qcddot1=zeros(size(qc));
tau1 = NewtonEuler(qc, qcdot1, qcddot1);


% centrifugal terms
tau2 = NewtonEuler(qc, qcdot, qcddot1);

% inertia terms
tau_inertial = tau- tau1-tau2;
tau3 = NewtonEuler(qc, qcdot1, qcddot);

figure;
hold on;
plot(time_span, tau(1,:), 'b');
plot(time_span, tau(2,:), 'r');
plot(time_span, tau(3,:), 'g');

title("compplete dynamics")
legend('Torque 1', 'Torque 2', 'Torque 3');
grid on
hold off

figure;
hold on;
plot(time_span, tau1(1,:), 'b');
plot(time_span, tau1(2,:), 'r');
plot(time_span, tau1(3,:), 'g');
title("gravity terms")
legend('Torque 1', 'Torque 2', 'Torque 3');
grid on
hold off

figure;
hold on;
plot(time_span, tau2(1,:), 'b');
plot(time_span, tau2(2,:), 'r');
plot(time_span, tau2(3,:), 'g');
title("coriolis term")
legend('Torque 1', 'Torque 2', 'Torque 3');
grid on
hold off


figure;
hold on;
plot(time_span, tau_inertial(1,:), 'b');
plot(time_span, tau_inertial(2,:), 'r');
plot(time_span, tau_inertial(3,:), 'g');
title("inertia terms")
legend('Torque 1', 'Torque 2', 'Torque 3');
grid on
hold off

figure;
subplot(3,1,1)
hold on
plot(time_span, tau(1,:), 'r', 'LineWidth',2);
plot(time_span, tau1(1,:), 'b');
plot(time_span, tau2(1,:), 'g');
plot(time_span, tau_inertial(1,:), 'y');
title("joint1")
legend('complete', 'gravity', 'coriolis', 'moments');
hold off

subplot(3,1,2)
hold on
plot(time_span, tau(2,:), 'r', 'LineWidth',2);
plot(time_span, tau1(2,:), 'b');
plot(time_span, tau2(2,:), 'g');
plot(time_span, tau_inertial(2,:), 'y');
grid on
legend('complete', 'gravity', 'coriolis', 'moments');
title("joint2")
hold off

subplot(3,1,3)
hold on
plot(time_span, tau(3,:), 'r', 'LineWidth',2);
plot(time_span, tau1(3,:), 'b');
plot(time_span, tau2(3,:), 'g');
plot(time_span, tau_inertial(3,:), 'y');
% legend('complete', 'gravity', 'coriolis', 'moments');
title("joint3")
hold off

function [Q]= NewtonEuler(qc, qcdot, qcddot)

xlabel('x')
ylabel('y')
zlabel('z')
view([0 0 1]);
axis equal;
xlim([-2 2]);
ylim([-2 2]);
zlim([-0.5 0.5]);
% legend('Link', 'Joint', 'Center of Mass')

% -------------------------------------------------------------------------
% the DH/ robot parameters
N_DOFS = 3;
dh.theta = [0 0 0];
dh.alpha = [-pi/2 pi/2 0];
dh.offset = [0 0 0];
dh.d = [0 1 1];
dh.a = [0 0 0];
dh.type = ['r' 'r' 'p'];

% -------------------------------------------------------------------------
% Rigid body paramaters: inertia, mass, and cener of ma\\ss
rb.I =  repmat(eye(3), 1, 1, N_DOFS);
rb.m = [1 1 1];


rb.r = [-0.25 0 0; -0.25 0 0; -0.25 0 0]';

EE = zeros(3, N_DOFS+1);
COM = zeros(3, N_DOFS);    

Q = invdyn(dh, rb, qc, qcdot, qcddot, [0; -9.8; 0]);



end


% this function calculates the inverse dynamics-------------------------------------------------------------------------
function Q = invdyn(dh, rb, qc, qcdot, qcddot, grav)
% Inverse dynamic with recursive Newton-Euler

if nargin < 6
    grav = [0;0;0];
end

z0 = [0; 0; 1];



for k = 1 : length(qc)  
    q = qc(k, :);
    qdot = qcdot(k, :);
    qddot = qcddot(k, :);

    N_DOFS = length(q);
    
    % ---------------------------------------------------------------------
    % Forward recursion
    for i = 1 : N_DOFS
        T = calc_transformation(i-1, i, dh, q);
        R = T(1:3, 1:3);
        p = [dh.a(i); dh.d(i)*sin(dh.alpha(i)); dh.d(i)*cos(dh.alpha(i))];
        
        if i > 1
            w(:, i) =  R'*(w(:, i-1) + z0.*qdot(i));
            wdot(:, i) = R'*(wdot(:, i-1) +  z0.*qddot(i) + ...
                cross(w(:, i-1), z0.*qdot(i)));
            vdot(:,i) = R'*vdot(:,i-1) + cross(wdot(:,i), p) + ...
                cross(w(:,i), cross(w(:,i),p));
        else
            w(:, i) =  R'*(z0.*qdot(i));
            wdot(:, i) = R'*(z0.*qddot(i));
            vdot(:,i) = R'*grav + cross(wdot(:,i), p) + ...
                cross(w(:,i), cross(w(:,i),p));
        end
    end
    
    % Dynamic simulation
    % Backward recursion starts here
    for i = N_DOFS:-1:1
        p = [dh.a(i); dh.d(i)*sin(dh.alpha(i)); dh.d(i)*cos(dh.alpha(i))];
        
        vcdot = vdot(:,i) + cross(wdot(:,i),rb.r(:,i)) + ...
            cross(w(:,i),cross(w(:,i),rb.r(:,i)));
        
        F = rb.m(i)*vcdot;
        N = rb.I(:,:,i)*wdot(:,i)+cross(w(:,i),rb.I(:,:,i)*w(:,1));
        
        if i < N_DOFS
            T = calc_transformation(i, i+1, dh, q);
            R = T(1:3, 1:3);
            n(:,i) = R*(n(:,i+1) + cross(R'*p, f(:,i+1))) + ...
                cross(rb.r(:,i)+p,F) + N;
            f(:,i) = R*f(:,i+1) + F;
        else
            n(:,i) = cross(rb.r(:,i)+p,F) + N;
            f(:,i) = F;
        end
        
        T = calc_transformation(i-1, i, dh, q);
        R = T(1:3, 1:3);
        
        if dh.type(i) == 'r'
            Q(i,k) = f(:,i)'*R'*z0;
        elseif dh.type(i) == 'p'
            Q(i,k) = n(:,i)'*R'*z0;
        end
    end
end
end

% -------------------------------------------------------------------------
function  T = calc_transformation(from, to, dh, qc)
% Transformation from one joint to another joint
% 0<=from<N_DOFS
% 0<to<=N_DOFS

T = eye(4);
N_DOFS = length(qc);

% check
if (from >= N_DOFS) || (from < 0) || (to <= 0) || (to >  N_DOFS)
    return;
end
% calculating forward kinematics to obtain rotational matrices
for i = from+1 : to
    if dh.type(i) == 'r'
        dh.theta(i) = qc(i);
    elseif dh.type(i) == 'p'
        dh.d(i) = qc(i);
    end
    
    ct = cos(dh.theta(i) + dh.offset(i));
    st = sin(dh.theta(i) + dh.offset(i));
    ca = cos(dh.alpha(i));
    sa = sin(dh.alpha(i));
    
    T = T * [ ct    -st*ca   st*sa     dh.a(i)*ct ; ...
              st    ct*ca    -ct*sa    dh.a(i)*st ; ...
              0     sa       ca        dh.d(i)    ; ...
              0     0        0         1          ];
end

end
