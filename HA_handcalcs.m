%% HA Aerothermal Calculations for Airframe Skin Temps %%
clear
clc

M_inf = [2:0.01:4.5];    % Freestream local mach Number (unitless)
Altitude = [10000:40:20000]; % Altitude of rocket (m)

[T_inf,~,P_inf,rho_inf] = atmosisa(Altitude);
P_abs =  P_inf ./101325; % Freestream pressure [atm] 

%k_titanium = 1078*10^3; % Thermal conductivity of Titanium (J/g*K)
rocket_length = 3;  % Length of rocket (m)
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

    for j = 1:1:length(T_local)
        k_pran(j) = k_air(round(T_local(j))); % Thermal Conductivity of air (Prandtl Number calculation) (W/ m^2 * K)
    end
    for k = 1:1:length(T_local)
        cp_pran(k) = cp(round(T_local(k))); % Specific heat of air at constant pressure (Prandtl Number calculation) (J/kg*K)
    end

    for y = 1:1:length(T_local)
        g_pran(y) = cp(round(T_local(y)))./(cp(round(T_local(y)))-R); % Ratio of Specfic Heats (unitless)
    end

for i = 1:1:length(M_inf)
    u_0 = 1.716*10^(-5);    % Reference dynamic viscosity (kg/m*s)
    T_pran_ref = 273.11;    % Sutherland reference temperature (K)
    S = 110.56;             % Sutherland's constant (unitless)

    u = u_0 * (T_local(i) / T_pran_ref) ^ (3 / 2) * (T_pran_ref + S) / (T_local(i) + S) * 10 ^ 3; % Dynamic viscosity of air (kg/m*s)

    Pr = u .* cp_pran(i) ./ k_pran(i);    % Prandtl Number (m^2/s)


    %% Adiabatic Wall Temperature
    r = Pr ^ (1 / 3);    % Turbulent Prandtl number (m^2/s)
    T_aw = T_local(i) * (1 + r * (g_pran(i) - 1) / 2 * m_cone(i) ^ 2);
    T_aw_stag = T_local(i) * (1 + (g_pran(i) - 1) / 2 * m_cone(i) ^ 2);

    %% Heat Transfer Coefficient with eador-Sfart Reference Temperature Method for Turbulent Flow

    T_star = T_local(i) * (0.55 * (1 + T_aw / T_local(i)) + 0.16 * r * ((g_pran(i) - 1) / 2) * m_cone(i) ^ 2);
    % T_star = T_inf + 0.5 * (T_w - T_inf) + 0.22 * (T_aw - T_w);

    for l = 1:1:length(T_local)
        k_star(l) = k_air(round(T_local(l))); % Thermal Conductivity of air (Star) (W/ m^2 * K)
    end

    for m = 1:1:length(T_local)
        cp_star(m) = cp(round(T_local(m))); % Specific heat of air at constant pressure (Star) (J/kg*K)
    end

    % Prandtl Number with Reference Temp
    u_star = u_0 * (T_star / T_pran_ref) ^ (3 / 2) * (T_pran_ref + S) / (T_star + S); % Dynamic viscosity of air (kg/m*s)

    Pr_star = u_star * cp_star(i) / k_star(i); % Reference temperature Prandtl number (m^2/s)

    % Reynolds Number
    R_diff = 0.028964917; % Molar mass of air (kg/mol)
    rho_star = P_local(i) / (R_diff * T_star);  % Reference temperature density (kg/m)
    u_star = u_0 * (T_star / T_pran_ref) ^ (3 / 2) * (T_pran_ref + S) / (T_star + S) * 10 ^ 3; % Dynamic viscosity of air (kg/m*s)

    Reynolds_star = rho_star * V_local(i) * rocket_length / u_star;

    % Nusselt Number Correlation
    Nu_star = 0.036 * Reynolds_star ^ 0.8 * Pr_star ^ (1 / 3);

    % HTC
    h(i) = k_star(i) * Nu_star / rocket_length; % Heat transfer coefficient (W/m^2*K)


    %% Heat Flux
    q(i) = h(i) * (T_aw - T_w);
end

%% HA Airframe Skin Temps %%

time = 1:1:251; % Time step
T_o = 300; % Initial Wall Temperature [K]
[t,y] = ode23s(@solver,time,T_o,odeset('RelTol',1e-12,'AbsTol',1e-12),h,T_inf,T_local,time);

hold on
grid on
plot(time,y); % Skin Temp Plot 
hold off

function Odesolver = solver(t,y,h,T_inf,T_local,time)
    for i = 1:1:length(time)
        rho = 4420; % Density of titanium [kg/m^3]
        V = 50; % Voulume of the airframe [m^3]
        A = 100; % Surface area of the airframe [m^2]
        c = 554.284; % Specific heat of titanium [J/kg * K]
        sigma = 5.6703e-8; % Stefan-Boltzmann Constant [W/m^2 * K^4]
        epi = 0.31; % Emmisivity
        dydt1 = (-A .* (h(i).*(y(1)-T_inf(i))+(sigma.*epi.*((y(1).^4) - (T_local(i)).^4)))) ./ (rho .* V .* c);
        Odesolver = dydt1;
    end
end
 


