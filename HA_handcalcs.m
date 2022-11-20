%% HA Aerothermal Calculations for Airframe Skin Temps %%
M_inf = [1];    % Freestream local mach Number (unitless)
Altitude = [1]; % Altitude of rocket (m)

[T_inf,~,P_inf,rho_inf] = atmosisa(Altitude);

k_titanium = 1078*10^3; % Thermal conductivity of aluminum (J/g*K)
rocket_length = 3;  % Length of rocket (m)
T_w = 300;      % Initial wall temp (K)

%% Gas Properties Lookup Tables 

% Tables for air cp, k as a function of temperature (Leons thing)

cp = lookup(T_inf);

%% Oblique Shock Relations

theta_cone = 30;    % Enter half-cone angle
g = 1.4;            % Enter the value of ratio of specific heats

%isentropic relations
To_T_freestream =1+(((g-1)/2).*M_inf.^2);
Po_P_freestream =(1+(((g-1)/2).*M_inf.^2)).^(g/(g-1));
rho0_rho_freestream =(1+(((g-1)/2).*M_inf.^2)).^(1/(g-1));

%freestream stagnation temperature and pressure(To = constant)
Po1 = Po_P_freestream*P_inf;
To = To_T_freestream*T_inf;

%calculation of shock wave angle(Q.10.1(a))
theta_shock = shock_angle(M_inf,theta_cone,g);

%calculates and returns the velocity and normal component of mach number
%along each ray between shock wave angle and the cone surface
[v,mn1]= flow_properties_behind_shock(M_inf,theta_shock,g);

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
rho_c = rho_freestream*((Pc/P_inf)^(1/g));

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
Po2_b=Po2_Po1_b*Po_P_behind_shock*P_inf;

%%pressure,tempearture and density immediately behind shock wave(Q.10(b))
Pb = Po2_b/Po_P_behind_shock;
Tb = To/To_T_behind_shock;
rho_b = rho_freestream*((Pb/P_inf)^(1/g));


T_local = 300;  % FreestreaM_inf local teM_infperature (K)
M_inf_local = 3;    % FreestreaM_inf local M_infach NuM_infber (unitless)
P_local = 101.325*10^6;  % FreestreaM_inf local absolute pressure (g/M_inf*s^2)
V_inf = 13; % FreestreaM_inf local velocity (M_inf/s)

%% Prandtl Number calculation
u_0 = 1.716*10^(-5);    % Reference dynamic viscosity (kg/M_inf*s)
T_pran_ref = 273.11;    % Sutherland reference temperature (273.11)
S = 110.56;             % Sutherland's constant (unitless)

cp = 1.005;  % Specific heat of air at constant pressure (J/g*K)
u = u_0 * (T_inf / T_pran_ref) ^ (3 / 2) * (T_pran_ref + S) / (T_inf + S) * 10 ^ 3; % Dynamic viscosity of air (g/M_inf*s)
k_air = 0.04441;  % Thermal conductivity of air (J/g*K)

Pr = u * cp / k_air;    % Prandtl Number (M_inf^2/s)


%% Adiabatic Wall TeM_infperature
r = Pr ^ (1 / 3);    % Turbulent Prandtl number (m^2/s)
T_aw = T_inf * (1 + r * (gamma - 1) / 2 * M_inf ^ 2);
T_aw_stag = T_inf * (1 + (gamma - 1) / 2 * M_inf ^ 2);

%% Heat Transfer Coefficient with eador-Sfart Reference Temperature Method for Turbulent Flow

T_star = T_inf * (0.55 * (1 + T_aw / T_inf) + 0.16 * r * ((gamma - 1) / 2) * M_inf ^ 2);
% T_star = T_inf + 0.5 * (T_w - T_inf) + 0.22 * (T_aw - T_w); 

% Prandtl Number with Reference Temp
cp_star = 1.216;  % Specific heat of air at constant pressure (J/g*K)
u_star = u_0 * (T_star / T_pran_ref) ^ (3 / 2) * (T_pran_ref + S) / (T_star + S); % Dynamic viscosity of air (g/M_inf*s)
k_air_star = 0.11007;  % Thermal conductivity of air (J/g*K)

Pr_star = u_star * cp_star / k_air_star; % Reference temperature Prandtl number (M_inf^2/s)

% Reynolds Number
R = 28.964917; % Specific gas constant (g/mol)
rho_star = P_inf / (R * T_star);  % Reference temperature density (g/m)
u_star = u_0 * (T_star / T_pran_ref) ^ (3 / 2) * (T_pran_ref + S) / (T_star + S) * 10 ^ 3; % Dynamic viscosity of air (g/M_inf*s)

Reynolds_star = rho_star * V_inf * rocket_length / u_star; 

% Nusselt Number Correlation
Nu_star = 0.036 * Reynolds_star ^ 0.8 * Pr_star ^ (1 / 3); 

% HTC
h = k_air_star * Nu_star / rocket_length; % Heat transfer coefficient (W/M_inf^2*K)


%% Heat Flux
q = h * (T_aw - T_w);


