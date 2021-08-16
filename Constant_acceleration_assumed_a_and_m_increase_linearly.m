%Create a meshgrid for the field
x1 = 15; 
y1 = 15;
z1 = 10730;
xa = linspace(0.1/1000,2.5/1000,x1);
ya = linspace(0.1/1000,2.5/1000,y1);
za = linspace(0.1/1000,1,z1);
[x,y,z]=meshgrid(xa,ya,za);
f = 750*(10^6);   %RF Frequency in Hz
V = 35*10^3;    %Vane Voltage in V
o = f/(2*pi);    %Angular Frequency 
lambda = (3*10^8)/f;    %Wavelength of RF in m
mass = 938.28*10^6;    %Mass of particles (protons) in beam in eV
initialenergy = 20*10^3;     %Initial kinetic energy of particles in eV
finalenergy = 2*10^6;    %Desired kinetic energy in eV
beta_i = sqrt(1-((mass/(initialenergy+mass))^2));    %Initial velocity of particles/speed of light
beta_f = sqrt(1-((mass/(finalenergy+mass))^2));     %Desired velocity of particles/speed of light
%Cell Length Values (Assumes a constant cell length increase)
L = zeros(x1,y1,z1);   
Li = beta_i*lambda;    %Length of modulation period (Based on velocity at start of cell)
Lf = beta_f*lambda;    %Length of last modulation period 
Lm = (Li+Lf)/2;    %Mean modulation period length
CN = 1/Lm;    %Number of modulation periods
deltaCL = (Lf-Li)/CN;    %modulation period length increase
CL = Li;    %CL = length of current modulation period
TL = CL;    %TL = distance along accelerator (z value)
%Modulation
m = zeros(x1,y1,z1);
mi = 1.7;    %Initial Modulation
mf = 2.5;    %Final Modulation
deltam = (mf-mi)/CN;    %Change in Modulation
mz = mi;
%Aperture 
a = zeros(x1,y1,z1);    
ai = 2/1000;    %Initial aperture size in m
af = 1/1000;    %Final aperture size in m
deltaa = (af-ai)/CN;    %Change in aperture size in m
az = ai;
%Loop cycles through every point in L,m and a matrices and puts new values
%in based on which 'modulation period' it's in.
    for l = 1:z1
        for j = 1:y1
            for i = 1:x1
                if (z(i,j,l) > TL-CL) && (z(i,j,l) <= TL)
                    L(i,j,l) = CL;
                    m(i,j,l) = mz;
                    a(i,j,l) = az;
                elseif (z(i,j,l) <= 1)
                    CL = CL + deltaCL;
                    mz = mz + deltam;
                    az = az + deltaa;
                    TL = TL + CL;
                    L(i,j,l) = CL;
                    m(i,j,l) = mz;
                    a(i,j,l) = az;
                else
                end
            end 
        end
    end
k = ((2*pi)/L);    
r = (abs(x.*x)+abs(y.*y)).^0.5;    %Radius from beam centre
G = ((besseli(0,k.*a))+(besseli(0,k.*m.*a)))./(((m.^2).*(besseli(0,k.*a)))+(besseli(0,m.*k.*a))); %Focusing efficiency (labelled as G instead of X to avoid confusion with x coordinate)
A = ((m.*m)-1)./(((m.*m).*(besseli(0,k.*a)))+(besseli(0,m.*k.*a)));    %Accelerating efficiency 
%Field components
Ex = -(((G*V)./(a.^2)).*x)-((k.*A*V)/2).*(besseli(1,k.*r)).*(x./r).*cos(k.*z);
Ey = (((G*V)./(a.^2)).*y)-((k.*A*V)/2).*(besseli(1,k.*r)).*(y./r).*cos(k.*z);
Ez = ((k.*A*V)/2).*(besseli(0,k.*r)).*sin(k.*z);
%clear f a V o x1 y1 z1 CL delta Li k r  A u v w q i j l  x y z beta CellNumber deltaCL lambda L
save('rf-track-2.0\Field_Map')    %Saves the workspace components