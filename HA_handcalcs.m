%% HA Aerothermal Calculations for Airframe & Maybe Nosecone Skin Temps %%
clear
clc

data = readmatrix('Important_stuff.xlsx');
Altitude_temp = data(:,1); % Altitude of rocket (m)
M_inf_temp = data(:,2); % Freestream local mach Number (unitless)
Max_mach_temp = max(M_inf_temp);
Max_mach_index = find(M_inf_temp == Max_mach_temp);

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

[T_inf,a_inf,P_inf,rho_inf] = atmosisa(Altitude);
P_abs =  P_inf ./101325; % Freestream pressure [atm] 

k_mat = 20.4; % Thermal conductivity of fiber glass E (W/m^2 * k)
rho_mat = 2563; % Density of fiber glass E (kg/m^3)
c_mat = 800; % Cp of fiber glass E (J / kg * K)
rocket_length = 1.778;  % Length of rocket (m)
T_w = 300;      % Initial wall temp (K)
R_nosecone = 0.04953; % Nosecone radius [m]
D_0 = 0.1016; % Outer Diameter of the Airframe [m]
D_i = 0.09525; % Inner Diameter of the Airframe [m]
V = (pi/4) * rocket_length*((D_0^2) - (D_i^2)); % Volume of the airframe [m^3]
A = pi * D_0 * rocket_length; % Surface area of the airframe [m^2]
Thermal_diff = k_mat / (rho_mat * c_mat); % Thermal Diffusivity of fiber glass E [m^2 / s]
Thick = (D_0 - D_i) / 2; % Thickness through which conduction occure through [m]
 

%% Gas Properties Lookup Tables 
R = 287.1; % Specific gas constant of air (J/ kg * K)

cp = airProp2(100:1:2500,'cp');
k_air = airProp2(100:1:2500,'k');
for i = 1:1:length(T_inf)
    g(i) = cp(round(T_inf(i)))./(cp(round(T_inf(i)))-R); % Ratio of Specfic Heats (unitless)
end

%% Oblique Shock Relations

theta_cone = 5.25;    % Half-cone angle
for i =1:1:length(M_inf)
    %isentropic relations
    To_T_freestream =1+(((g(i)-1)./2).*M_inf(i).^2);
    Po_P_freestream =(1+(((g(i)-1)./2).*M_inf(i).^2)).^(g(i)./(g(i)-1));
    rho0_rho_freestream =(1+(((g(i)-1)./2).*M_inf(i).^2)).^(1./(g(i)-1));

    %freestream stagnation temperature and pressure(To = constant)
    Po1 = Po_P_freestream.*P_abs(i);
    To(i) = To_T_freestream.*T_inf(i);

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
    m_cone=mach(length(mach));

    %isentropic relations
    To_T_surface=1+(((g(i)-1)./2).*m_cone.^2);
    Po_P_surface=(1+(((g(i)-1)./2).*m_cone.^2)).^(g(i)./(g(i)-1));
    rho0_rho_surface=(1+(((g(i)-1)./2).*m_cone.^2)).^(1./(g(i)-1));

    %ratio of stagnation pressure across the shock wave using equation()
    Po2_Po1=((((g(i)+1)./2.*mn1.^2)./(1+(((g(i)-1)./2).*mn1.^2))).^(g(i)./(g(i)-1)))./((((2.*g(i)./(g(i)+1)).*mn1.^2)-((g(i)-1)./(g(i)+1))).^(1/(g(i)-1)));
    Po2=Po2_Po1.*Po1;

    %pressure,tempearture and density on cone surface(Q.10(b))
    Pc = Po2./Po_P_surface;
    Tc = To(i)./To_T_surface;
    rho_c = rho_inf(i).*((Pc./P_abs(i)).^(1./g(i)));

    %the mach number immediately behind the shock wave will be the first
    %element in the array of mach number calculated above
    m_behind_shock(i) = mach(1);

    %isentropic relations
    To_T_behind_shock(i)=1+(((g(i)-1)./2).*m_behind_shock(i).^2);
    Po_P_behind_shock(i)=(1+(((g(i)-1)./2).*m_behind_shock(i).^2)).^(g(i)./(g(i)-1));
    rho0_rho_behind_shock(i)=(1+(((g(i)-1)./2).*m_behind_shock(i).^2)).^(1./(g(i)-1));

    %shock relations
    mn_b = m_behind_shock(i).*sin(theta_shock*pi./180);
    Po2_Po1_b=((((g(i)+1)./2.*mn_b.^2)./(1+(((g(i)-1)./2).*mn_b.^2))).^(g(i)./(g(i)-1)))./((((2.*g(i)./(g(i)+1)).*mn_b.^2)-((g(i)-1)/(g(i)+1))).^(1/(g(i)-1)));
    Po2_b=Po2_Po1_b.*Po_P_behind_shock(i).*P_abs(i);

    %%pressure,tempearture and density immediately behind shock wave(Q.10(b))
    P_local(i) = (Po2_b./Po_P_behind_shock(i)) * 101325;
    T_local(i) = To(i)./To_T_behind_shock(i); % Freestream local temperature (K)
    rho_local(i) = rho_inf(i).*((P_local(i)./P_abs(i)).^(1./g(i)));
    V_local(i) = m_behind_shock(i) * sqrt(R*T_local(i)*g(i)); % Freestream local velocity (m/s)
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
    T_aw(i) = T_local(i) * (1 + r * (g_pran - 1) / 2 * m_behind_shock(i) ^ 2);
    T_aw_stag = T_local(i) * (1 + (g_pran - 1) / 2 * m_behind_shock(i) ^ 2);

    %% Heat Transfer Coefficient with eador-Sfart Reference Temperature Method for Turbulent Flow

    T_star = T_local(i) * (0.55 * (1 + T_aw(i) / T_local(i)) + 0.16 * r * ((g_pran - 1) / 2) * m_behind_shock(i) ^ 2);
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

%% Sutton Graves
k = 1.7415e-4;
for i=1:1:length(M_inf)
    V_inf = M_inf(i) * a_inf(i);
    q_s(i) = k * sqrt(rho_inf(i)/R_nosecone) * (V_inf ^ 3);
end

%% HA Airframe Skin Temps %%
t = 1:1:200; % Time step
biot = (h / k_mat) * (V / A);
Fourier = (Thermal_diff * t) / (Thick^2);
biot_2 = (h / k_mat) * (Thick);

IC = 300;

Biot_known = [0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100];
C_known = [1.0025 1.005 1.0075 1.0099 1.0124 1.0148 1.0173 1.0197 1.0222 1.0246 1.0365 1.0483 1.0598 1.0712 1.0932 1.1143 1.1345 1.1539 1.1724 1.1902 1.2071 1.3384 1.4191 1.4698 1.5029 1.5253 1.5411 1.5526 1.5611 1.5677 1.5919 1.5973 1.5993 1.6002 1.6015];
Zeta_known = [0.1412 0.1995 0.2440 0.2814 0.3143 0.3438 0.3709 0.3960 0.4195 0.4417 0.5376 0.6170 0.6856 0.7465 0.8516 0.9408 1.0184 1.0873 1.1490 1.2048 1.2558 1.5994 1.7887 1.9081 1.9898 2.0490 2.0937 2.1286 2.1566 2.1795 2.2881 2.3261 2.3455 2.3572 2.3809];
C = interp1(Biot_known,C_known,biot_2);
zeta = interp1(Biot_known,Zeta_known,biot_2);


for i = 1:1:length(t)
    if biot < 0.1 % Lum sum
        [t,T_internal] = ode45(@(t,T_internal) energy_equation(t, T_internal, h(i), T_aw(i), V, A, T_local(i),rho_mat,c_mat), t, IC);
    elseif Fourier(i) > 0.2 % 1st term approx
        T_internal(i) = (IC - T_local(i)) * (C(i) * exp((-zeta(i)^2)*Fourier(i))) + T_local(i);
    elseif Fourier(i) <= 0.2
        Part1 = erfc(Thick/(2*sqrt(Thermal_diff*t(i))));
        Part2 = exp((h(i)*Thick/k_mat) + ((h(i)^2) * Thermal_diff * t(i))/(k_mat^2));
        Part3 = erfc((Thick/(2*sqrt(Thermal_diff*t(i)))) + (h(i)*sqrt(Thermal_diff*t(i)) / k_mat));
        T_internal(i) = (Part1 - (Part2 * Part3))*(T_local(i) - IC) + IC;
        %T_skin(i) = ((2*q(i)*sqrt((Thermal_diff*time(i))/pi)/k_mat)*exp(-(Thick^2)/(4*Thermal_diff*time(i)))) - ((q(i)*Thick/k_mat) * (erfc(Thick/(2*sqrt(Thermal_diff*time(i)))))) + IC;  
        %T_skin(i) = IC;
        %T_skin(i) =(T_local(i) - IC) * (erfc(Thick/(2*Thermal_diff*time(i))) - exp(((h(i)*Thick)/k_mat) + (h(i)^2) * Thermal_diff * time(i)/(k_mat^2)) * erfc((Thick/(2*sqrt(Thermal_diff*time(i)))) + (h(i)*sqrt(Thermal_diff*time(i))/k_mat))) + IC;
    end
end

%IC = 300;
%[time_sol,T_skin] = ode45(@(t,T_skin) energy_equation(t, T_skin, h, T_aw, time, V, A), time, IC);


yyaxis right
plot(t, T_internal, 'linewidth', 2.5)
title('Performance')
xlabel("Time (s)")
ylabel('Internal Temperature (K)')
yyaxis left
plot(t, M_inf, 'linewidth', 2.5)
ylabel('Mach Number')
grid on;

% plot thermal quantities
figure

sgtitle('Thermal Variables') 

subplot(3,3,1)
plot(t, h, 'linewidth', 2.5)
title('Heat Transfer Coeff (Airframe)')
xlabel('Time (s)')
ylabel('h (W/m^2K)')
grid on;

subplot(3,3,2)
plot(t, q, 'linewidth', 2.5)
title('Heat Flux (Airframe)')
xlabel('Time (s)')
ylabel('q (W/m^2)')
grid on;

%subplot(3,2,3)
%plot(time, T_aw, 'linewidth', 2.5)
%title('Adiabatic Wall Temperature')
%xlabel('Time (s)')
%ylabel('T_a_w (K)')
%grid on;

subplot(3,3,3)
plot(t, biot, 'linewidth', 2.5)
title('Biot Number (Airframe with Fiberglass)')
xlabel('Time (s)')
ylabel('Bi')
grid on;

subplot(3,3,4)
plot(t,q_s,'linewidth',2.5)
title('Heat flux of Stagnation Point')
xlabel('Time (s)')
ylabel('q (W/m^2)')
grid on;

subplot(3,3,5)
plot(t,Fourier,'linewidth',2.5)
title('Fourier Number (Fiberglass)')
xlabel('Time (s)')
ylabel('F_o')
grid on;

subplot(3,3,6)
plot(t, biot_2, 'linewidth', 2.5)
title('Biot Number_2 (Airframe with Fiberglass)')
xlabel('Time (s)')
ylabel('Bi')
grid on;

subplot(3,3,7)
plot(t,T_inf-T_internal, 'linewidth', 2.5)
title('Delta Temp (T_i_n_f - T_i_n_t_e_r_n_a_l) (Airframe)')
xlabel('Time (s)')
ylabel('Delta Temp (K)')
grid on; 


function T_sdot = energy_equation(time_sol, T_internal, h, T_aw, V, A, T_local,rho_mat,c_mat)
    sigma = 5.6703e-8; % Stefan-Boltzmann Constant [W/m^2 * K^4]
    epi = 0.31; % Emmisivity
    %h = interp1(t,h,time_sol);
    %T_aw = interp1(t,T_aw,time_sol);
    %T_local = interp1(t,T_local,time_sol);

    T_sdot = A * (h * (T_aw - T_internal) - sigma * epi * (T_internal ^ 4 - T_local ^ 4)) / (rho_mat * V * c_mat);
end





