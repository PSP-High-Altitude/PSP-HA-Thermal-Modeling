% The following code provides a solution to solve the Taylor Maccoll
% equation for supersonic flow over a cone. 
% The user needs to specify the freestream Mach number(M), the half cone angle(theta_cone) and the value
% of ratio of specific heat(g). 
% The solver specified below calls 4 sub-routines to solve the entire
% flowfield. 

% By: Rohan Kokate
% Reference: Fundamentals of Aerodynamics, J. D. Anderson.

clear;
clc;

M = 5;              %enter freestream mach no
theta_cone = 50;    %enter half-cone angle
g = 1.4;            %enter the value of ratio of specific heats

%flow properties upstream of shock(standard sea level conditions)
P_freestream = 1;
T_freestream = 287;
rho_freestream = 1.225;

%isentropic relations
To_T_freestream =1+(((g-1)/2).*M.^2);
Po_P_freestream =(1+(((g-1)/2).*M.^2)).^(g/(g-1));
rho0_rho_freestream =(1+(((g-1)/2).*M.^2)).^(1/(g-1));

%freestream stagnation temperature and pressure(To = constant)
Po1 = Po_P_freestream*P_freestream;
To = To_T_freestream*T_freestream;

%calculation of shock wave angle(Q.10.1(a))
theta_shock = shock_angle(M,theta_cone,g);

%calculates and returns the velocity and normal component of mach number
%along each ray between shock wave angle and the cone surface
[v,mn1]= flow_properties_behind_shock(M,theta_shock,g);

%total velocity value along each ray
velocity=sqrt((v(:,1).^2)+(v(:,2).^2));

%mach number based on equation()
mach=sqrt(2./(((velocity.^(-2))-1).*(g-1)));

%mach number on cone surface is the last value in the array of mach numbers
%calculated above
m_cone=mach(length(mach));

%isentropic relations
To_T_surface=1+(((g-1)/2).*m_cone.^2);
Po_P_surface=(1+(((g-1)/2).*m_cone.^2)).^(g/(g-1));
rho0_rho_surface=(1+(((g-1)/2).*m_cone.^2)).^(1/(g-1));

%ratio of stagnation pressure across the shock wave using equation()
Po2_Po1=((((g+1)/2.*mn1.^2)./(1+(((g-1)/2).*mn1.^2))).^(g/(g-1)))./((((2*g/(g+1)).*mn1.^2)-((g-1)/(g+1))).^(1/(g-1)));
Po2=Po2_Po1*Po1;

%pressure,tempearture and density on cone surface(Q.10(b))
Pc = Po2/Po_P_surface;
Tc = To/To_T_surface;
rho_c = rho_freestream*((Pc/P_freestream)^(1/g));

%the mach number immediately behind the shock wave will be the first
%element in the array of mach number calculated above
m_behind_shock = mach(1);

%isentropic relations
To_T_behind_shock=1+(((g-1)/2).*m_behind_shock.^2);
Po_P_behind_shock=(1+(((g-1)/2).*m_behind_shock.^2)).^(g/(g-1));
rho0_rho_behind_shock=(1+(((g-1)/2).*m_behind_shock.^2)).^(1/(g-1));

%shock relations
mn_b = m_behind_shock*sin(theta_shock*pi/180);
Po2_Po1_b=((((g+1)/2.*mn_b.^2)./(1+(((g-1)/2).*mn_b.^2))).^(g/(g-1)))./((((2*g/(g+1)).*mn_b.^2)-((g-1)/(g+1))).^(1/(g-1)));
Po2_b=Po2_Po1_b*Po_P_behind_shock*P_freestream;

%%pressure,tempearture and density immediately behind shock wave(Q.10(b))
Pb = Po2_b/Po_P_behind_shock;
Tb = To/To_T_behind_shock;
rho_b = rho_freestream*((Pb/P_freestream)^(1/g));

%Calculation of drag co-efficient based on equation()(Q.10.3)
C_D = (2*((Pc/P_freestream)-1))/(g*M^2);

% End of routine