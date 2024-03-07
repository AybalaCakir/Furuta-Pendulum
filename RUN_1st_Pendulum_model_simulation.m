Jm = 3.87e-7;
J = 2.982e-5;
Jb = 0.0044;
kg = 70;
km = 0.00767;
Rm = 2.6;
m = 0.14;
l = 0.165;
r = 0.175;
g = 9.81;

J1 = Jm*kg^2 + J;

tfin = 1;
Y = 1;
t = (0:0.001:tfin)';


A = [0, 0, 1, 0; 0, 0, 0, 1; g*(Jb+m*r^2+J1)/(l*(J1+Jb)), 0, 0, r*km^2*kg^2/(l*Rm*(J1+Jb)); -m*r*g/(J1+Jb), 0, 0, -km^2*kg^2/(Rm*(J1+Jb))];
B = [0, 0, -r*km*kg/(l*Rm*(J1+Jb)), km*kg/(Rm*(J1+Jb))]';
C_orig = [0 1 0 0];
C = [l r 0 0];

% u = -(1/2)*km*kg*Y*pi^3*cos(pi*t/tf)/tf^3-(1/2)*Y*(Rm*l*pi^4*J1+Rm*l*pi^4*Jb-g^2*Rm*tf^4*Jb-g^2*Rm*tf^4*m*r^2-g^2*Rm*tf^4*J1)*sin(pi*t/tf)/(r*km*kg*g*tf^4)-(1/2)*Y*(-g^2*Rm*tf^4*m*r^2-g^2*Rm*tf^4*Jb-g^2*Rm*tf^4*J1)/(r*km*kg*g*tf^4)

% u = -(1/2)*Y*(Rm*l*sin((1/2)*pi*(2*t-tf)/tf)*pi^4*J1+Rm*l*sin((1/2)*pi*(2*t-tf)/tf)*pi^4*Jb-g^2*Rm*tf^4*Jb-g^2*Rm*tf^4*Jb*sin((1/2)*pi*(2*t-tf)/tf)-g^2*Rm*tf^4*m*r^2-g^2*Rm*tf^4*m*r^2*sin((1/2)*pi*(2*t-tf)/tf)-g^2*Rm*tf^4*J1-g^2*Rm*tf^4*J1*sin((1/2)*pi*(2*t-tf)/tf)+cos((1/2)*pi*(2*t-tf)/tf)*pi^3*tf*km^2*kg^2*r*g)/(r*km*kg*g*tf^4);

% u = -6*Y*km*kg*t.^2/(r*tf^3)+(-(12*(g^2*Rm*Jb+g^2*Rm*m*r^2+g^2*Rm*J1))*Y/(r*km*kg*g^2*tf^3)+6*Y*km*kg/(r*tf^2))*t+(6*(g^2*Rm*Jb+g^2*Rm*m*r^2+g^2*Rm*J1))*Y/(r*km*kg*g^2*tf^2)+12*l*Y*km*kg/(g*r*tf^3)

u = -20*Y*km*kg*t.^4/(r*tfin^5)+(-(80*(g^2*Rm*Jb+g^2*Rm*m*r^2+g^2*Rm*J1))*Y/(r*km*kg*g^2*tfin^5)+60*Y*km*kg/(r*tfin^4))*t.^3+((180*(g^2*Rm*Jb+g^2*Rm*m*r^2+g^2*Rm*J1))*Y/(r*km*kg*g^2*tfin^4)+(240*l*Y/(g*r*tfin^5)-60*Y/(r*tfin^3))*km*kg)*t.^2+(-(480*(-Rm*l*J1-Rm*l*Jb))*Y/(r*km*kg*g*tfin^5)-(120*(g^2*Rm*Jb+g^2*Rm*m*r^2+g^2*Rm*J1))*Y/(r*km*kg*g^2*tfin^3)+(-360*l*Y/(g*r*tfin^4)+20*Y/(r*tfin^2))*km*kg).*t+(360*(-Rm*l*J1-Rm*l*Jb))*Y/(r*km*kg*g*tfin^4)+(20*(g^2*Rm*Jb+g^2*Rm*m*r^2+g^2*Rm*J1))*Y/(r*km*kg*g^2*tfin^2)+120*l*Y*km*kg/(g*r*tfin^3);

sys = ss(A,B,C_orig,0);

lsim(sys,u,t);
ysoll = Y*(10*t.^2/tfin^2-20*t.^3/tfin^3+15*t.^4/tfin^4-4*t.^5/tfin^5);

plot(t,u); grid on;
plot(t,ysoll);

thetasoll = Y*(10*t.^2/tfin^2-20*t.^3/tfin^3+15*t.^4/tfin^4-4*t.^5/tfin^5)/r-l*Y*(20/tfin^2-120.*t/tfin^3+180*t.^2/tfin^4-80*t.^3/tfin^5)/(g*r);
plot(t,thetasoll);
