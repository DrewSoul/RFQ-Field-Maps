%Create meshgrid for field
x1 = 15;
y1 = 15;
z1 = 10730;
xa = linspace(0.1/1000,2.5/1000,15);
ya = linspace(0.1/1000,2.5/1000,15);
za = linspace(0.1/1000,1,10730);
[x,y,z]=meshgrid(xa,ya,za);
f = 750*(10^6);   %RF Frequency in Hz
V = 35*10^3;    %Vane Voltage in V
o = f/(2*pi);    %Angular Frequency
lambda = (3*10^8)/f;    %Wavelength of RF in m
mass = 938.28*10^6;    %Mass of particles in beam (Proton) in eV
InitialEnergy = 20*10^3;    %Initial energy of particle beam in eV
FinalEnergy = 2*10^6; %Desired energy of particle beam in eV
Energy = InitialEnergy;
CN = 1102.185855;    %Number of modulation periods
deltaE = (FinalEnergy - InitialEnergy)/CN;    %(Constant) energy increaese per modulation period
CL = lambda*sqrt(1-((mass/(InitialEnergy+mass))^2));    %modulation period length (Initial) 
TL = CL;    %Length along accelerator (z coordinate)
%Modulation
m = zeros(x1,y1,z1);
mi = 1.7;    %Initial Modulation
mf = 2.5;    %Final Modulation
deltam = (mf-mi)/CN;    %Change in Modulation
mz = mi;
%Aperture 
a = zeros(x1,y1,z1);    
ai = 2/1000;    %Initial aperture size in m
af = 0.5/1000;    %Final aperture size in m
deltaa = (af-ai)/CN;    %Change in aperture size
az = ai;
%Loop cycles through every point in L,m and a matrices and puts new values
%in based on which modulation period it's in.
   
    for l = 1:z1
        for j = 1:y1
            for i = 1:x1
                if (z(i,j,l) > TL-CL) && (z(i,j,l) <= TL)
                    L(i,j,l) = CL;
                    m(i,j,l) = mz;
                    a(i,j,l) = az;
                elseif (z(i,j,l) <= 1)
                    Energy = Energy + deltaE;
                    CL = lambda*sqrt(1-((mass/(Energy+mass))^2));
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