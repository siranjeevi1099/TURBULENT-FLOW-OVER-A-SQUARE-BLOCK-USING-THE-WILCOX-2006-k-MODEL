function [] = bound()
% Purpose: Specify boundary conditions for a calculation

% constants
global NPI NPJ U_IN YMAX Cmu Ti low_bound high_bound 
% variables
global y u v T m_in m_out y_v F_u k eps f cx cy x y LARGE SP omega mu rho omegaw mueff

% Set velocity at the inlet


% step_width  = 0.05;   % m
% step_height = 0.01;   % m
% 
% for I = 1:NPI+1
% 
%     if I == 1
%         % leftmost cell (no I-1 available) → use face coordinate directly
%         cx = x(1);
%     else
%         % interior cells → average of two faces
%         cx = 0.5*(x(I) + x(I-1));
%     end
% 
%     for J = 1:NPJ+1
%         if J == 1
%             % bottommost cell (no J-1 available)
%             cy = y(1);
%         else
%             cy = 0.5*(y(J) + y(J-1));
%         end
% 
%         % check if inside the step block
%         if (cx <= step_width) && (cy <= step_height)
%             u(I,J)=0;   % solid region
%         end
%     end
% end

u(2,1:NPJ+1) = U_IN;

% u(2,ceil((NPJ+2)/3):NPJ+1) = U_IN;
% u(2,1:ceil((NPJ+2)/3))=0;

% Set k and eps at the inlet


% Set k and eps at the inlet
k(1,3:NPJ)     = 1.5*(U_IN*Ti)^2; % at inlet
k(:,1:2) = 0;     
k(:,NPJ+1:NPJ+2) = 0;
omegaw(1,3:NPJ)   = Cmu^0.75 *k(1,3:NPJ).^0.5/(0.07*YMAX*0.5); % at inlet
omegaw(:,1:2) = 60* (2.E-5/1.2) / (0.075 * (YMAX/NPJ)) ;
omegaw(:,NPJ+1:NPJ+2) = 60* (2.E-5/1.2) / (0.075 * (YMAX/NPJ)) ;

%for step
% k(1,ceil((NPJ+2)/3):NPJ+2)     = 1.5*(U_IN*Ti)^2; % at inlet
% k(1,1:ceil((NPJ+2)/3))     = 1e-3; % at STEP
% eps(1,ceil((NPJ+2)/3):NPJ+2)   = Cmu^0.75 *k(1,ceil((NPJ+2)/3):NPJ+2).^1.5/(0.07*YMAX*0.5); % at inlet
% eps(1,1:ceil((NPJ+2)/3))     = 1e-4; % at STEP

% Fix temperature at the walls in Kelvin
T(1:NPI+2,1) = 273.; % bottom wall
T(1:NPI+2,NPJ+2) = 273.; % top wall
   


% begin: globcont();=======================================================
% Purpose: Calculate mass in and out of the calculation domain to correct for the continuity at outlet.
convect();

m_in = 0.;
m_out = 0.;

for J = 2:NPJ+1
    j = J;
    AREAw = y_v(j+1) - y_v(j); % See fig. 6.3
    m_in  = m_in  + F_u(2,J)*AREAw;
    m_out = m_out + F_u(NPI+1,J)*AREAw;
end
% end: globcont()==========================================================

% Velocity and temperature gradient at outlet = zero:
% Correction factor m_in/m_out is used to satisfy global continuity
u(NPI+2,2:NPJ+1) = u(NPI+1,2:NPJ+1);%m_in/m_out;
v(NPI+2,2:NPJ+1) = v(NPI+1,2:NPJ+1);
k(NPI+2,2:NPJ+1) = k(NPI+1,2:NPJ+1);
eps(NPI+2,2:NPJ+1) = eps(NPI+1,2:NPJ+1);
T(NPI+2,1:NPJ+2) = T(NPI+1,1:NPJ+2);
end
