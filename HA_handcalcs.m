%% HA Aerothermal Calculations for Airframe Skin Temps %%
clear
clc

M_inf = [2 3 5];    % Freestream local mach Number (unitless)
Altitude = [1000 5000 10000]; % Altitude of rocket (m)

[T_inf,~,P_inf,rho_inf] = atmosisa(Altitude);
P_abs =  P_inf ./101300; 

k_titanium = 1078*10^3; % Thermal conductivity of Titanium (J/g*K)
rocket_length = 3;  % Length of rocket (m)
T_w = 300;      % Initial wall temp (K)

%% Gas Properties Lookup Tables 
R = 0.028964917; % Specific gas constant (kg/mol)
% Tables for air cp, k as a function of temperature (Leons thing)

cp = 1010; %Cp (J/kg * K)
g = [1.1 1.2 1.3]; % Ratio of Specfic Heats (unitless)
%% Oblique Shock Relations

theta_cone = 40;    % Half-cone angle
for i =1:1:length(M_inf)
    %isentropic relations
    To_T_freestream =1+(((g(i)-1)./2).*M_inf(i).^2);
    Po_P_freestream =(1+(((g(i)-1)./2).*M_inf(i).^2)).^(g(i)./(g(i)-1));
    rho0_rho_freestream =(1+(((g(i)-1)./2).*M_inf(i).^2)).^(1./(g(i)-1));

    %freestream stagnation temperature and pressure(To = constant)
    Po1 = Po_P_freestream(i).*P_abs(i);
    To = To_T_freestream(i).*T_inf(i);

    %calculation of shock wave angle(Q.10.1(a))
    theta_shock = shock_angle(M_inf(i),theta_cone,g(i));

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
    To_T_surface=1+(((g-1)./2).*m_cone.^2);
    Po_P_surface=(1+(((g-1)./2).*m_cone.^2)).^(g./(g-1));
    rho0_rho_surface=(1+(((g-1)./2).*m_cone.^2)).^(1./(g-1));

    %ratio of stagnation pressure across the shock wave using equation()
    Po2_Po1=((((g+1)./2.*mn1.^2)./(1+(((g-1)./2).*mn1.^2))).^(g./(g-1)))./((((2.*g./(g+1)).*mn1.^2)-((g-1)./(g+1))).^(1/(g-1)));
    Po2=Po2_Po1.*Po1;

    %pressure,tempearture and density on cone surface(Q.10(b))
    Pc = Po2./Po_P_surface;
    Tc = To./To_T_surface;
    rho_c = rho_inf.*((Pc./P_abs).^(1./g));

    %the mach number immediately behind the shock wave will be the first
    %element in the array of mach number calculated above
    m_behind_shock = mach(1);

    %isentropic relations
    To_T_behind_shock=1+(((g-1)./2).*m_behind_shock.^2);
    Po_P_behind_shock=(1+(((g-1)./2).*m_behind_shock.^2)).^(g./(g-1));
    rho0_rho_behind_shock=(1+(((g-1)./2).*m_behind_shock.^2)).^(1./(g-1));

    %shock relations
    mn_b = m_behind_shock.*sin(theta_shock*pi./180);
    Po2_Po1_b=((((g+1)./2.*mn_b.^2)./(1+(((g-1)./2).*mn_b.^2))).^(g./(g-1)))./((((2.*g./(g+1)).*mn_b.^2)-((g-1)/(g+1))).^(1/(g-1)));
    Po2_b=Po2_Po1_b.*Po_P_behind_shock.*P_abs;

    %%pressure,tempearture and density immediately behind shock wave(Q.10(b))
    P_local = Po2_b./Po_P_behind_shock;
    T_local = To./To_T_behind_shock; % FreestreaM_inf local temperature (K)
    rho_local = rho_inf.*((P_local./P_abs).^(1./g));
end

M_inf_local = 3;    % Freestream local mach Number (unitless)
V_inf = 13; % Freestream local velocity (m/s)
%% Prandtl Number calculation
u_0 = 1.716*10^(-5);    % Reference dynamic viscosity (kg/m*s)
T_pran_ref = 273.11;    % Sutherland reference temperature (K)
S = 110.56;             % Sutherland's constant (unitless)

u = u_0 * (T_inf / T_pran_ref) ^ (3 / 2) * (T_pran_ref + S) / (T_inf + S) * 10 ^ 3; % Dynamic viscosity of air (kg/m*s)
k_air = 5;  % Thermal conductivity of air (J/kg*K)

Pr = u * cp / k_air;    % Prandtl Number (m^2/s)


%% Adiabatic Wall Temperature
r = Pr ^ (1 / 3);    % Turbulent Prandtl number (m^2/s)
T_aw = T_inf * (1 + r * (g - 1) / 2 * M_inf ^ 2);
T_aw_stag = T_inf * (1 + (g - 1) / 2 * M_inf ^ 2);

%% Heat Transfer Coefficient with eador-Sfart Reference Temperature Method for Turbulent Flow

T_star = T_inf * (0.55 * (1 + T_aw / T_inf) + 0.16 * r * ((g - 1) / 2) * M_inf ^ 2);
% T_star = T_inf + 0.5 * (T_w - T_inf) + 0.22 * (T_aw - T_w); 

% Prandtl Number with Reference Temp
cp_star = 5;  % Specific heat of air at constant pressure (J/kg*K)
u_star = u_0 * (T_star / T_pran_ref) ^ (3 / 2) * (T_pran_ref + S) / (T_star + S); % Dynamic viscosity of air (kg/m*s)
k_air_star = 5;  % Thermal conductivity of air (J/kg*K)

Pr_star = u_star * cp_star / k_air_star; % Reference temperature Prandtl number (m^2/s)

% Reynolds Number
rho_star = P_inf / (R * T_star);  % Reference temperature density (kg/m)
u_star = u_0 * (T_star / T_pran_ref) ^ (3 / 2) * (T_pran_ref + S) / (T_star + S) * 10 ^ 3; % Dynamic viscosity of air (kg/m*s)

Reynolds_star = rho_star * V_inf * rocket_length / u_star; 

% Nusselt Number Correlation
Nu_star = 0.036 * Reynolds_star ^ 0.8 * Pr_star ^ (1 / 3); 

% HTC
h = k_air_star * Nu_star / rocket_length; % Heat transfer coefficient (W/m^2*K)


%% Heat Flux
q = h * (T_aw - T_w);


