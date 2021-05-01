%% Calculating eigen frequencies of a ring resonator with densiry defect and
%% defect of midline radius. Cross section of resonator - rectangle.
%% Bernoulli beam FE implemented.

%% Preparation
clc
clear all
close all

%% Initial data
N_e = 100;                  % Number of finite elements in model
R = 20e-3;                  % Radius of midline of ring resonator (nominal value)
b = 3e-3;                   % Height of rasonator cross section
h = 0.5e-3;                 % Thickness of ring resonator (nominal value)
E = 7e10;                   % Young modulus (nominal value)
ro = 2210;                  % Density  (nominal value)

%% Defects initial data. Each defect distributed by two harmonics.
% Amplitudes of harmonics
dR = 0*2e-3*[0,1];          % Midline radius defect
db = 0*[1,1];               % Height defect
dh = 0*0.25e-3*[1,0];       % Thickness defect
dE = 0*[1,1];               % Young modulus defect
dro = [1,1]*ro/100;       % Density defect
% Numbers of harmonics
kR = [1,3];
kb = [1,3];
kh = [1,3];
kE = [1,3];
kro = [8,12];
% Initial angle of harmonics [degrees]
fi_R = [0,0]*pi/180;
fi_b = [0,0]*pi/180;
fi_h = [0,0]*pi/180;
fi_E = [0,0]*pi/180;
fi_ro = [0,0]*pi/180;

scl = R/10;                 % Printing scale for eigen vectors
scl_dim = 1e3;              % Dimension printing scale

% Arrays with defects for FE mesh
R_err = zeros(N_e+1,1);
b_err = zeros(N_e,1);
h_err = zeros(N_e,1);
E_err = zeros(N_e,1);
ro_err = zeros(N_e,1);

%% Calculation
N_n = N_e+1;                % Number of nodes
fi = (0:2*pi/N_e:2*pi)';    % Polar angle coordinates of nodes
fi_mid = fi(1:end-1) + (fi(2)-fi(1))/2; % Polar angle coordinates of midpoints of elements

% Defects calculation
for num=1:length(dR)
    R_err = R_err + dR(num)*cos(kR(num)*(fi-fi_R(num)));
end
for num=1:length(db)
    b_err = b_err + db(num)*cos(kb(num)*(fi_mid-fi_b(num)));
end
for num=1:length(dh)
    h_err = h_err + dh(num)*cos(kh(num)*(fi_mid-fi_h(num)));
end
for num=1:length(dE)
    E_err = E_err + dE(num)*cos(kE(num)*(fi_mid-fi_E(num)));
end
for num=1:length(dro)
    ro_err = ro_err + dro(num)*cos(kro(num)*(fi_mid-fi_ro(num)));
end

% Sin and Cos of angle fi
S = sin(fi); C = cos(fi);
% Coordinates of perfect ring resonator (used for printing)
x0 = R*C;
y0 = R*S;
% Coordinates of imperfect coordinates of ring resonator FE nodes
% (used for calculation stiffness and mass matrices)
x = R*C + R_err.*C;
y = R*S + R_err.*S;

%% Stiffness and mass matrices forming
K_mat = zeros(3*N_e);
M_mat = zeros(3*N_e);

for i=1:N_e
    if (i~=N_e) 
        % if onsidered element is not the last one
        El_num = i;     % element number
        n1 = El_num;    % 1-st node number of element
        n2 = El_num+1;  % 2-nd node number of element
        Dof_num = (3*(n1-1)+1:3*(n1-1)+6);  % indeces of DOFs present in considering element
        % E_t - is a structure with information about element
        E_t.x1 = x(n1);
        E_t.y1 = y(n1);
        E_t.x2 = x(n2);
        E_t.y2 = y(n2);
        E_t.b = b+b_err(El_num);
        E_t.h = h+h_err(El_num);
        E_t.E = E+E_err(El_num);
        E_t.p = ro+ro_err(El_num);
        
        % Calculation of stiffness and mass matrices of an element
        [K_e, M_e] = Str_Beam_KM(E_t);
        % Assembling matrices of considered element to the whole resonator
        % matrices
        K_mat(Dof_num,Dof_num) = K_mat(Dof_num,Dof_num) + K_e;
        M_mat(Dof_num,Dof_num) = M_mat(Dof_num,Dof_num) + M_e;
    else
        % if considered element is the last one it connects to the first
        % finite element as resonator is closed in circumferential direction
        El_num = i;
        n1 = El_num;
        n2 = 1;         % for the last element 2-nd node is the 1-st node of the mesh
        Dof_num = [(3*(i-1)+1:3*(i-1)+3), (3*(n2-1)+1:3*(n2-1)+3)];
        
        E_t.x1 = x(n1); % see above
        E_t.y1 = y(n1);
        E_t.x2 = x(n2);
        E_t.y2 = y(n2);
        E_t.b = b+b_err(El_num);
        E_t.h = h+h_err(El_num);
        E_t.E = E+E_err(El_num);
        E_t.p = ro+ro_err(El_num);
        
        [K_e, M_e] = Str_Beam_KM(E_t);
        K_mat(Dof_num,Dof_num) = K_mat(Dof_num,Dof_num) + K_e;
        M_mat(Dof_num,Dof_num) = M_mat(Dof_num,Dof_num) + M_e;
    end
end

% Calculating of eigenvalues
[u1,p] = eig(K_mat,M_mat);
% Getting sqrt so that p is circular frequency now
p = sqrt(diag(p));
% Sorting frequencies
[p, ind] = sort(p);
u1 = u1(:,ind);
% Excluding zero frequencies corresponding to movements as a solid boby
p = p(4:end);
u1 = u1(:,4:end);

% Getting frequencies corresponding to wineglass forms
% (elliptical modes, modes with k=2 circumferential wavenumber)
k = 2;
% Theoretical value of not perturbated frequency
p_th = sqrt(E*h^2/(12*ro*R^4)*k.^2.*(k.^2-1).^2./(k.^2+1));

col = ['g';'r';'m';'k'];

N = 1;

% Printing geometry of a resonator
figure(1)
axis equal;
grid on;
hold on;
plot(scl_dim*x0,scl_dim*y0,'k--','linewidth',1,'Marker','none');
plot(scl_dim*x,scl_dim*y,'b','linewidth',2,'Marker','s');
legend('resonator FE model')
xlabel('x, [mm]','FontName','TimesNewRoman','fontsize',36)
ylabel('y, [mm]','FontName','TimesNewRoman','fontsize',36)
set(gca,'FontName','TimesNewRoman','fontsize',36)
plot(xlim,[0 0],'k');
plot([0 0],ylim,'k');

%% Getting and plotting eigen modes
for num = 1:1:2
    
    i = 1:N_e;          % indeces of DOFs
    iu = 3*(i-1)+1;     % indeces of u displacements
    iv = 3*(i-1)+2;     % indeces of v displacements
    itz = 3*(i-1)+3;    % indeces of rotational angle theta
    
    u = u1(iu,N+num-1); u = [u; u(1)]; % u displacement
    v = u1(iv,N+num-1); v = [v; v(1)]; % v displacement
    tz = u1(itz,N+num-1); tz = [tz; tz(1)]; % rotational angle theta
    
    % Scaling eigen modes
    delta_R = sqrt(u.^2+v.^2);  % Radial diasplacement 
    [n, m] = max(abs(delta_R)); % Maximal radial displacement
    n = sign(delta_R(m))*n;     % Scale of naturel modes to unit
    u = scl*u/n;                % scaling u disp
    v = scl*v/n;                % scaling v disp
    
    % New coordinates of ring resonator for mode with number "num"
    X = x+u;
    Y = y+v;
    % Plotting eigen modes
    figure(num+1)
    axis equal;
    grid on; hold on;
    plot(scl_dim*X,scl_dim*Y,col(num),'linewidth',2,'Marker','none')
    
    xlabel('x, [mm]','FontName','TimesNewRoman','fontsize',36)
    ylabel('y, [mm]','FontName','TimesNewRoman','fontsize',36)
    set(gca,'FontName','TimesNewRoman','fontsize',36)
    xlim = get(gca, 'XLim');
    ylim = get(gca, 'YLim');
    plot(xlim,[0 0],'k');
    plot([0 0],ylim,'k');
    plot(scl_dim*x,scl_dim*y,'b--','linewidth',1);
    
end
figure(2)
legend('"heavy" mode')
figure(3)
legend('"light" mode')


dp = abs(p(2)-p(1));                        % Numerical FEM frequency split
dp_th = 3*dro(1)*dro(2)/10/2*p_th/ro^2;     % Theoretical frequency split

nn = min(kro);
Fn = (516*nn + 929*nn^2 + 720*nn^3 + 290*nn^4 + 60*nn^5 + 5*nn^6)/(90 + 36*nn);
dp_th = (2*(3+2*nn)/5/Fn + 0.2*(nn==1)*(abs(diff(kro))~=4))/2*p_th*dro(1)*dro(2)/ro^2;

display(['Numerical frequency split, Hz     ',num2str(dp/2/pi)])
display(['Theoretical frequency split, Hz   ',num2str(dp_th/2/pi)])







