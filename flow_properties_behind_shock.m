
%exactly the same function as cone_angle.m except that this function
%returns the value of velocity and normal component of mach number for each
%ray behind the shock wave. This basically solves the entire flow field
%behind the shock wave from the shock wave to the cone surface. The reason
%for not using function cone_angle.m itslef is that, the velocity and mach
%number values are useless in the bisectional method calculation of shock
%wave

function [v,mn1]= flow_properties_behind_shock(m1,theta_shock,g)
global gamma
g=gamma;
theta_shock=theta_shock.*(pi)/180;
delta=atan(2.*cot(theta_shock).*(((m1.^2).*(sin(theta_shock).^2)-1)./((m1.^2).*(gamma+cos(2*theta_shock))+2)));
mn1=m1.*sin(theta_shock);
mn2=sqrt(((mn1.^2)+(2/(gamma-1)))./((2*gamma./(gamma-1)).*(mn1.^2)-1));
m2=mn2./sin(theta_shock-delta);
v_in=(2/((gamma-1).*(m2.^2))+1).^(-.5);
v_rin=v_in.*cos(theta_shock-delta);
v_thin=v_in.*sin(theta_shock-delta);
endtheta=0.1.*(pi)/180;
options=odeset('Events','on');
[theta,v]=ode23('taylor_maccoll_eqn',[theta_shock, endtheta],[v_rin, v_thin], options);
theta_degree=theta.*(180/(pi));
cone_angle =theta_degree(length(theta_degree));

end
% End of subroutine