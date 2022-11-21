% Import look up table instructions as individual column vectors
% Create Simulink Model using instructions

% "CpLookup" should be renamed to Simulink Model file name
load_system("CpLookup");

%simin should be replaced with whatever independent var
simin = 2200;
sim("CpLookup", simin);
out = ans.simout;
disp(out);
