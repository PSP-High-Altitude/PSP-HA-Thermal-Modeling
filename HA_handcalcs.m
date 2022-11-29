%% HA Aerothermal Calculations for Airframe & Maybe Nosecone Skin Temps %%
clear
clc

data = readmatrix('Important_stuff.xlsx');
Altitude_temp = data(:,1); % Altitude of rocket (m)
M_inf_temp = data(:,2); % Freestream local mach Number (unitless)
Max_mach_temp = max(M_inf_temp);
Max_mach_index = find(M_inf_temp == Max_mach_temp);
% Apogee = max(Altitude_temp);
% Apogee_index = find(Altitude_temp == Apogee);

Start_mach = find((M_inf_temp >= 2));
Start_mach_index = min(Start_mach);
count = 1;

for i = Start_mach_index:1:Max_mach_index
        M_inf_temp2(count) = M_inf_temp(i);
        Altitude_temp2(count) = Altitude_temp(i);
        count = count + 1;
end

indexVal = find(M_inf_temp2 == max(M_inf_temp2));
count = 1;

for i = round(linspace(1,indexVal,200))
   M_inf(count) = M_inf_temp2(i);
   Altitude(count) = Altitude_temp2(i);
   count = count + 1;
end

% M_inf = linspace(3,3,200);   % Freestream local mach Number (unitless)
% Altitude = linspace(25908,25908,200); % Altitude of rocket (m)

[T_inf,~,P_inf,rho_inf] = atmosisa(Altitude);
P_abs =  P_inf ./101325; % Freestream pressure [atm] 

k_titanium = 20.4; % Thermal conductivity of Titanium (W/m^2 * k)
rocket_length = 1.88;  % Length of rocket (m)
T_w = 300;      % Initial wall temp (K)

%% Gas Properties Lookup Tables 
R = 287.1; % Specific gas constant of air (J/ kg * K)
% Tables for air cp, k as a function of temperature (Leons thing)

cp = airProp2(100:1:2500,'cp');
k_air = airProp2(100:1:2500,'k');
for i = 1:1:length(T_inf)
    g(i) = cp(round(T_inf(i)))./(cp(round(T_inf(i)))-R); % Ratio of Specfic Heats (unitless)
end

%% Oblique Shock Relations

theta_cone = 40;    % Half-cone angle
for i =1:1:length(M_inf)
    %isentropic relations
    To_T_freestream =1+(((g(i)-1)./2).*M_inf(i).^2);
    Po_P_freestream =(1+(((g(i)-1)./2).*M_inf(i).^2)).^(g(i)./(g(i)-1));
    rho0_rho_freestream =(1+(((g(i)-1)./2).*M_inf(i).^2)).^(1./(g(i)-1));

    %freestream stagnation temperature and pressure(To = constant)
    Po1 = Po_P_freestream.*P_abs(i);
    To = To_T_freestream.*T_inf(i);

    %calculation of shock wave angle(Q.10.1(a))
    theta_shock = shock_angle(M_inf(i),theta_cone,g(i));

    %calculates and returns the velocity and normal component of mach number
    %along each ray between shock wave angle and the cone surface
    [v,mn1]= flow_properties_behind_shock(M_inf(i),theta_shock,g(i));

    %total velocity value along each ray
    velocity=sqrt((v(:,1).^2)+(v(:,2).^2));

    %mach number based on equation()
    mach=sqrt(2./(((velocity.^(-2))-1).*(g(i)-1)));

    %mach number on cone surface is the last value in the array of mach numbers
    %calculated above
    m_cone(i)=mach(length(mach));

    %isentropic relations
    To_T_surface=1+(((g(i)-1)./2).*m_cone(i).^2);
    Po_P_surface=(1+(((g(i)-1)./2).*m_cone(i).^2)).^(g(i)./(g(i)-1));
    rho0_rho_surface=(1+(((g(i)-1)./2).*m_cone(i).^2)).^(1./(g(i)-1));

    %ratio of stagnation pressure across the shock wave using equation()
    Po2_Po1=((((g(i)+1)./2.*mn1.^2)./(1+(((g(i)-1)./2).*mn1.^2))).^(g(i)./(g(i)-1)))./((((2.*g(i)./(g(i)+1)).*mn1.^2)-((g(i)-1)./(g(i)+1))).^(1/(g(i)-1)));
    Po2=Po2_Po1.*Po1;

    %pressure,tempearture and density on cone surface(Q.10(b))
    Pc = Po2./Po_P_surface;
    Tc = To./To_T_surface;
    rho_c = rho_inf(i).*((Pc./P_abs(i)).^(1./g(i)));

    %the mach number immediately behind the shock wave will be the first
    %element in the array of mach number calculated above
    m_behind_shock = mach(1);

    %isentropic relations
    To_T_behind_shock=1+(((g(i)-1)./2).*m_behind_shock.^2);
    Po_P_behind_shock=(1+(((g(i)-1)./2).*m_behind_shock.^2)).^(g(i)./(g(i)-1));
    rho0_rho_behind_shock=(1+(((g(i)-1)./2).*m_behind_shock.^2)).^(1./(g(i)-1));

    %shock relations
    mn_b = m_behind_shock.*sin(theta_shock*pi./180);
    Po2_Po1_b=((((g(i)+1)./2.*mn_b.^2)./(1+(((g(i)-1)./2).*mn_b.^2))).^(g(i)./(g(i)-1)))./((((2.*g(i)./(g(i)+1)).*mn_b.^2)-((g(i)-1)/(g(i)+1))).^(1/(g(i)-1)));
    Po2_b=Po2_Po1_b.*Po_P_behind_shock.*P_abs(i);

    %%pressure,tempearture and density immediately behind shock wave(Q.10(b))
    P_local(i) = (Po2_b./Po_P_behind_shock) * 101325;
    T_local(i) = To./To_T_behind_shock; % Freestream local temperature (K)
    rho_local(i) = rho_inf(i).*((P_local(i)./P_abs(i)).^(1./g(i)));
    V_local(i) = m_cone(i) * sqrt(R*T_local(i)*g(i)); % Freestream local velocity (m/s)
end

%% Prandtl Number calculation
for i = 1:1:length(M_inf)
    
    k_pran = k_air(round(T_local(i))); % Thermal Conductivity of air (Prandtl Number calculation) (W/ m^2 * K)
 
    cp_pran = cp(round(T_local(i))); % Specific heat of air at constant pressure (Prandtl Number calculation) (J/kg*K)

    g_pran = cp(round(T_local(i)))./(cp(round(T_local(i)))-R); % Ratio of Specfic Heats (unitless)

    u_0 = 1.716*10^(-5);    % Reference dynamic viscosity (kg/m*s)
    T_pran_ref = 273.11;    % Sutherland reference temperature (K)
    S = 110.56;             % Sutherland's constant (unitless)

    u = u_0 * (T_local(i) / T_pran_ref) ^ (3 / 2) * (T_pran_ref + S) / (T_local(i) + S); % Dynamic viscosity of air (kg/m*s)

    Pr = u .* cp_pran ./ k_pran;    % Prandtl Number (m^2/s)


    %% Adiabatic Wall Temperature
    r = Pr ^ (1 / 3);    % Turbulent Prandtl number (m^2/s)
    T_aw(i) = T_local(i) * (1 + r * (g_pran - 1) / 2 * m_cone(i) ^ 2);
    T_aw_stag = T_local(i) * (1 + (g_pran - 1) / 2 * m_cone(i) ^ 2);

    %% Heat Transfer Coefficient with eador-Sfart Reference Temperature Method for Turbulent Flow

    T_star = T_local(i) * (0.55 * (1 + T_aw(i) / T_local(i)) + 0.16 * r * ((g_pran - 1) / 2) * m_cone(i) ^ 2);
    % T_star = T_inf + 0.5 * (T_w - T_inf) + 0.22 * (T_aw - T_w);

    k_star = k_air(round(T_star)); % Thermal Conductivity of air (Star) (W/ m^2 * K)

    cp_star = cp(round(T_star)); % Specific heat of air at constant pressure (Star) (J/kg*K)

    % Prandtl Number with Reference Temp
    u_star = u_0 * (T_star / T_pran_ref) ^ (3 / 2) * (T_pran_ref + S) / (T_star + S); % Dynamic viscosity of air (kg/m*s)

    Pr_star = u_star * cp_star / k_star; % Reference temperature Prandtl number (m^2/s)

    % Reynolds Number
    rho_star = P_local(i) / (R * T_star);  % Reference temperature density (kg/m)
    u_star = u_0 * (T_star / T_pran_ref) ^ (3 / 2) * (T_pran_ref + S) / (T_star + S); % Dynamic viscosity of air (kg/m*s)

    Reynolds_star = rho_star * V_local(i) * rocket_length / u_star;

    % Nusselt Number Correlation
    Nu_star = 0.036 * Reynolds_star ^ 0.8 * Pr_star ^ (1 / 3);

    % HTC
    h(i) = k_star * Nu_star / rocket_length; % Heat transfer coefficient (W/m^2*K)


    %% Heat Flux
    q(i) = h(i) * (T_aw(i) - T_w);
end

%% HA Airframe Skin Temps %%

V = 0.0010586853; % Volume of the airframe [m^3]
A = 1.412706672; % Surface area of the airframe [m^2]
biot = h / k_titanium * V / A;

IC = 300;
time = 1:1:200; % Time step
[time_sol,T_skin] = ode45(@(t,T_skin) energy_equation(t, T_skin, h, T_aw, time, V, A), time, IC);


% plot skin temperature
yyaxis right
plot(time_sol, T_skin, 'linewidth', 2.5)
title('Performance')
xlabel("Time (s)")
ylabel('Skin Temperature (K)')
yyaxis left
plot(time, M_inf, 'linewidth', 2.5)
ylabel('Mach Number')
grid on;

% plot thermal quantities
figure

sgtitle('Thermal Variables') 

subplot(2,2,1)
plot(time, h, 'linewidth', 2.5)
title('Heat Transfer Coeff')
xlabel('Time (s)')
ylabel('h (W/m^2K)')
grid on;

subplot(2,2,2)
plot(time, q, 'linewidth', 2.5)
title('Heat Flux')
xlabel('Time (s)')
ylabel('q (W/m^2)')
grid on;

subplot(2,2,3)
plot(time, T_aw, 'linewidth', 2.5)
title('Adiabatic Wall Temperature')
xlabel('Time (s)')
ylabel('T_a_w (K)')
grid on;

subplot(2,2,4)
plot(time, biot, 'linewidth', 2.5)
title('Biot Number')
xlabel('Time (s)')
ylabel('Bi')
grid on;

function T_sdot = energy_equation(t, T_skin, h, T_aw, time, V, A)
    sigma = 5.6703e-8; % Stefan-Boltzmann Constant [W/m^2 * K^4]
    epi = 0.31; % Emmisivity
    h = interp1(time,h,t);
    T_aw = interp1(time,T_aw,t);
    rho = -0.1416*T_aw + 4461; % Density of titanium [kg/m^3]
    c = 0.1093*T_aw + 551.54; % Specific heat of titanium [J/kg * K]

    T_sdot = A * (h * (T_aw - T_skin) - sigma * epi * (T_skin ^ 4 - T_aw ^ 4)) / (rho * V * c);
end





