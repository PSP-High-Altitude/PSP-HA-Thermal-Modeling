%this function actually solves the state variable form of the taylor
%maccoll equation derived in section 2. The option routine set in the
%cone_angle.m and flow_properties_behind_shock.m is executed here. 


function [value, isterminal, direction] = taylor_maccoll_eqn(t,g,flag)
global gamma

gdot(1)=g(2);

a=(gamma-1)/2;

Numerator = a.*(2*((g(1)).^3-g(1))+(g(2)-g(2).*g(1).^2-g(2).^3).*cot(t)-2*g(1).*g(2).^2)-(g(1).^2.*g(2).^2);
Denominator = a.*((g(2)).^2+(g(1)).^2 - 1) + (g(2)).^2;
gdot(2) = Numerator./Denominator;
% this event detects when g(2)(angular velocity) is zero, stop the iteration from proceeding

if (nargin<3)||isempty(flag)

value=[gdot(1);gdot(2)];

%used youtube video(reference 3) to understand the functioning of ode23
else flag=='events';
value=g;            % values to check
isterminal=[0;1];   % terminal on g(2)
direction=[0;-1];   % detects falling slope(not sure what this does)
end

%end of subroutine