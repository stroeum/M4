function [vcn] = Collision_Frequency(c,n,nc,nn,Tc,Tn,h)

% Formulas from Schunk and Nagy, 2000, pp. 97--99, Tables 4.4, 4.5, 4.6 %
Tr = (Tc+Tn)/2.0;
if     (strcmp(c,'O2+')  && strcmp(n,'CO2')) 
    vcn = 5.63e-10*nn*1e-6;
elseif (strcmp(c,'CO2+') && strcmp(n,'CO2'))
    vcn = (Tr>850)*2.85e-11*(nn*1e-6)*sqrt(Tr)*(1-0.083*log10(Tr))^2.0;
elseif (strcmp(c,'O+')   && strcmp(n,'CO2')) 
    vcn = 8.95e-10*nn*1e-6;
elseif (strcmp(c,'e')    && strcmp(n,'CO2')) 
    vcn = 3.68e-8*(nn*1e-6)*(1+4.1e-11*abs(4500 - Tc)^2.93);
elseif (strcmp(c,'O2+')  && strcmp(n,'O'))   
    vcn = 2.31e-10*nn*1e-6;
elseif (strcmp(c,'CO2+') && strcmp(n,'O'))   
    vcn = 1.76e-10*nn*1e-6;
elseif (strcmp(c,'O+')   && strcmp(n,'O'))   
    vcn = (Tr>235)*3.67e-11*(nn*1e-6)*sqrt(Tr)*(1-0.064*log10(Tr))^2.0;
elseif (strcmp(c,'e')    && strcmp(n,'e'))   
    vcn = 54.6/sqrt(2)*(nn*1e-6)*Tc^-1.5;
elseif (strcmp(c,'e')    && strcmp(n,'O'))   
    vcn = 8.9e-11*(nn*1e-6)*(1+5.7e-4*Tc)*sqrt(Tc);
elseif (strcmp(c,'e')    && strcmp(n,'O2+'))   
    vcn = 54.6*(nn*1e-6)*Tc^-1.5;
elseif (strcmp(c,'e')    && strcmp(n,'CO2+'))
    vcn = 54.6*(nn*1e-6)*Tc^-1.5;
elseif (strcmp(c,'e')    && strcmp(n,'O+'))
    vcn = 54.6*(nn*1e-6)*Tc^-1.5;
elseif (strcmp(c,'O2+')  && strcmp(n,'O2+'))
    vcn = 0.16*(nn*1e-6)*Tn^-1.5;
elseif (strcmp(c,'CO2+') && strcmp(n,'O2+'))
    vcn = 0.12*(nn*1e-6)*Tn^-1.5;
elseif (strcmp(c,'O+')   && strcmp(n,'O2+'))
    vcn = 0.26*(nn*1e-6)*Tn^-1.5;
elseif (strcmp(c,'O2+')  && strcmp(n,'CO2+'))
    vcn = 0.17*(nn*1e-6)*Tn^-1.5;
elseif (strcmp(c,'O2+')  && strcmp(n,'O+'))
    vcn = 0.13*(nn*1e-6)*Tn^-1.5;
elseif (strcmp(c,'CO2+') && strcmp(n,'CO2+'))
    vcn = 0.14*(nn*1e-6)*Tn^-1.5;
elseif (strcmp(c,'CO2+') && strcmp(n,'O+'))
    vcn = 0.10*(nn*1e-6)*Tn^-1.5;
elseif (strcmp(c,'O+')   && strcmp(n,'CO2+'))
    vcn = 0.27*(nn*1e-6)*Tn^-1.5;
elseif (strcmp(c,'O+')   && strcmp(n,'O+'))
    vcn = 0.22*(nn*1e-6)*Tn^-1.5;
else   printf('Warning: unknown collision\n\tAssuming null collision frequency.\n'); vcn=0.0;
end
if (isnan(vcn) || nc<=0)
    sprintf('nu: error detected: h=%f, c=%d n=%d,Tc=%f,Tn=%f\n',h,c,n,Tc,Tn);
end
end

% PetscReal vcn(PetscInt c, PetscInt n, PetscReal h)
% {
%   PetscReal ne,nc,nn,Tc,Tn;
%
%   ne = Profile(7,h);
%
%   Tc = Profile(9,h);
%   if (strcmp(c,'O2+' || strcmp(c,'CO2+' || strcmp(c,'O+') {
%     if (strcmp(c,'O2+')  nc = 0.71*ne;
%     if (strcmp(c,'CO2+') nc = 0.25*ne;
%     if (strcmp(c,'O+')   nc = 0.04*ne;
%     Tc = Profile(9,h);
%   } else if (strcmp(c,'e') {
%     nc = ne;
%     Tc = Profile(8,h);
%   } else {
%     printf('Error: Unknown charged species\n');
%     exit(0);
%   }
%
%   if (strcmp(n,'CO2')  nn = Profile(5,h);
%   if (strcmp(n,'O')    nn = Profile(4,h);
%   Tn = Profile(10,h);
%
%   return Vcn(c,n,nc,nn,Tc,Tn,h);
% }