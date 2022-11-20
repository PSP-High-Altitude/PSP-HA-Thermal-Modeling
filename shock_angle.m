

%this sub-routine iteratures through a series of shock angles until it
%reaches the shock angle which is supported by the cone angle specified by
%the user
function theta_shock = shock_angle(M,theta_cone,g)


iterations_max = 60;    % put the maximum number of iterations
iteration = 0;          %initialise 

% There are two shock waves for a specific cone angle(strong and weak
% shock wave). Using equation(), solve for the strong shock angle.(2D
% wedge case)
theta_shock_max = 50;   
cone_angle_max = cone_angle(M,theta_shock_max,g);   %calculates cone angle which supports the strong shock wave

%calculate the weak shock wave, which is the mach wave(2D wedge case)
theta_shock_min = asin(1/M)*180/pi+0.2;
cone_angle_min = cone_angle(M,theta_shock_min,g);   %calculates cone angle which supports the weak shock wave


%setup iteration to check for the correct shock angle at the surface of the
%cone using the bisectional iteration method discussed in section 2

while (iteration<iterations_max)
theta_shock_mean = (theta_shock_max+theta_shock_min)/2;
cone_angle_mean = cone_angle(M,theta_shock_mean,g);
if cone_angle_mean>theta_cone
theta_shock_max=theta_shock_mean;
else
theta_shock_min=theta_shock_mean;
end
iteration = iteration+1;
end

%return the value of the correct shock wave angle
theta_shock=theta_shock_mean;


%End of subroutine