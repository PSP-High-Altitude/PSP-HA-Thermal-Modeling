

%this sub-routine solves for the cone angle for a specific shock wave
%angle. This function returns the value of cone angle which is then checked
%with the user defined value of cone angle in the previous
%function(shock_angle.m)

function cone_angle = cone_angle(M,theta_shock,g)

global gamma
gamma=g; 
m1=M;

% starting shock angle
theta_shock=theta_shock.*(pi)/180;

%here we use the oblique shock relations for calculations
% deflection angle just behind the wave
delta=atan(2.*cot(theta_shock).*(((m1.^2).*(sin(theta_shock).^2)-1)./((m1.^2).*(gamma+cos(2*theta_shock))+2)));

% normal component of the free stream Mach No.
mn1=m1.*sin(theta_shock);

% normal component of the Mach no. after the shock
mn2=sqrt(((mn1.^2)+(2/(gamma-1)))./((2*gamma./(gamma-1)).*(mn1.^2)-1));

% calculation of Mach no. after the shock
m2=mn2./sin(theta_shock-delta);

% initial dimensionless component of velocity after shock
v_initial=(2/((gamma-1).*(m2.^2))+1).^(-.5);

% radial and normal components of the velocity calculated
% as boundary conditions
v_r_initial = v_initial.*cos(theta_shock-delta);
v_t_initial = v_initial.*sin(theta_shock-delta);

% lower bound for cone angle (limiting case)(cannot put to be abs 0)
endtheta=0.1.*(pi)/180;     

% options command switches on event detection which is shown later(acts
% like a stopper when the tangential component of velocity becomes zero)
options=odeset('Events','on');

% we use ode23 solver which returns the numerical solution for the current
% shock angle(Runga-Kutta method)

[theta,v]=ode23('taylor_maccoll_eqn',[theta_shock, endtheta],[v_r_initial, v_t_initial], options); 

% converts the angle values to degrees
theta_degree=theta.*(180/(pi));

% returns the cone angle from the numerical values found
cone_angle =theta_degree(length(theta_degree));


% End of subroutine