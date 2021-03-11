clear all
%%

% Initialising


% Set up for arteries is as follows. ArtertyName = [1.Length, 2.AreaIn, 3.AreaOut, 4.Daughters Check, 5.DaughterNode1 or Total Resistance, 6.DaughterNode2 or Total Compliance, 7.Compl0, 8.Compl1, 9.Wall Thickness, 10. Youngs' Modulus, 11.Number of
% System]

AscendingAorta = [4E-2, 1.470E-2, 1.440E-2, 1, 2, 3, -0.4853E-06, 3.0794E-09, 0.164E-2, 4E+5 1]; %

AorticArchA = [2.0E-2, 1.120E-2, 1.120E-2, 1, 14, 15, 1.16650E-06, 2.8208E-09, 0.132E-2, 4E+5,2]; %

Innominate = [3.4E-2, 0.620E-2, 0.620E-2, 1, 4, 5, 4.9882E-06, 2.162E-09, 0.086E-2, 4E+5, 3]; %

RSubclavianA = [3.4E-2, 0.423E-2, 0.423E-2, 1, 6, 7, 7.15050E-6,1.1710E-09, 0.067E-2, 4E+5, 4]; %

RCartoid = [17.7E-2, 0.37E-2, 0.37E-2, 1, 12, 13, 7.74630E-6, 1.5746E-09, 0.063E-2, 4E+5, 5]; %

RVertebral = [14.8E-2, 0.188E-2,0.183E-2, 0, 0.60100E+10, 0.30955E-10, 7.6606E-06, 2.2096E-10, 0.046E-2, 8E+5,6];

RSubclavianB = [42.2E-2, 0.403E-2, 0.236E-2, 1, 8, 9, 9.2673E-6, 1.0976E-09, 0.066E-2, 4E+5, 7];

RRadial = [23.5E-2, 0.174E-2,0.142E-2, 0, 0.52800E+10, 0.3523E-10, 7.459E-6, 2.0325E-10, 0.044E-2, 8E+5, 8];

RUlnarA = [6.7E-2, 0.215E-2, 0.215E-2, 1, 10, 11, 8.0504E-6, 2.603E-10, 0.049E-2, 8E+5, 9];

RInterosseus = [7.9E-2, 0.091E-2, 0.091E-2, 0, 0.84300E+11, 0.22068E-11, 3.8843E-6, 3.8352E-11, 0.028E-2, 16E+5, 10];

RUlnarB = [17.1E-2, 0.203E-2,0.183E-2, 0, 0.52800E+10, 0.22069E-10, 7.8845E-6, 2.4264E-10, 0.047E-2, 8E+5, 11];

RInternalCartoid = [17.7E-2, 0.177E-2, 0.083E-2, 0, 0.13900E+11, 0.13384E-10, 6.8426E-6, 1.5751E-10, 0.045E-2, 8E+5, 12];

RExternalCartoid = [17.7E-2, 0.177E-2, 0.083E-2, 0, 0.13900E+11, 0.13384E-10, 6.8426E-6, 1.5751E-10, 0.045E-2, 8E+5, 13];

AorticArchB = [2.0E-2, 1.120E-2 1.120E-2, 1, 18, 19, 1.5509E-6, 2.7589E-09, 0.132E-2, 4E+5, 14];

LCartoid = [20.8E-2, 0.37E-2,0.37E-2, 1, 16, 17, 7.74680E-6, 1.5745E-09, 0.067E-2, 4E+5, 15];

LInternalCartoid = [17.7E-2, 0.177E-2,0.083E-2, 0, 0.13900E+11, 0.13384E-10, 6.8426E-6, 1.551E-10, 0.045E-2, 8E+5, 16];

LExternalCartoid = [17.7E-2, 0.177E-2,0.083E-2, 0, 0.13900E+11, 0.13384E-10, 6.8426E-6, 1.551E-10, 0.045E-2, 8E+5, 17];

ThoraicAortaA = [5.2E-2, 0.999E-2,0.999E-2, 1, 26, 27, 2.0236E-6, 2.6817E-09, 0.12E-2, 4E+5, 18];

LSubclavianA = [3.4E-2, 0.423E-2,0.423E-2, 1, 20, 21, 7.1505E-6, 1.717E-09, 0.067E-2, 4E+5, 19];

Vertebral = [14.8E-2, 0.188E-2,0.183E-2 , 0, 0.60100E+10, 0.30955E-10, 7.6606E-6, 2.2096E-10, 0.046E-2, 8E+5, 20];

LSubclavianB = [42.2E-2, 0.403E-2,0.236E-2, 1, 22, 23, 9.2673E-6, 1.0976E-09, 0.066E-2, 4E+5, 21];

LRadial = [23.5E-2, 0.174E-2,0.142E-2, 0, 0.52800E+10, 0.35235E-10, 7.459E-6, 2.0325E-10, 0.044E-2, 8E+5, 22];

LUlnarA = [6.7E-2, 0.215E-2,0.215E-2, 1, 24, 25, 8.0504E-6, 2.603E-10, 0.049E-2, 8E+5, 23];

LInterosseous = [7.9E-2, 0.091E-2,0.091E-2, 0, 0.84300E+11, 0.22068E-11, 3.8843E-6, 3.8352E-11, 0.028E-2, 16E+5, 24];

LUlnarB = [17.1E-2, 0.203E-2,0.1833E-2, 0, 0.52800E+10, 0.22069E-10, 7.8845E-6, 2.4264E-10, 0.047E-2, 8E+5, 25];

Intercostals = [8.0E-2, 0.2E-2,0.5E-2, 0, 0.13900E+10, 0.13384E-09, 5.543E-7, 2.918E-9, 0.05E-2, 4E+5, 26]; 

ThoraicAortaB = [10.4E-2, 0.675E-2, 0.645E-2, 1, 28, 29, 4.5981E-6, 2.2348E-09, 0.09E-2, 4E+5, 27];

AbdominalAortaA = [5.3E-2, 0.61E-2,0.61E-2, 1, 34, 35, 4.9553E-6, 2.1683E-09, 0.084E-2, 4E+5, 28];

CeliacA = [1E-2, 0.39E-2,0.39E-2, 1, 30, 31, 7.5619E-6, 1.6201E-09, 0.064E-2, 4E+5, 29];

CeliacB = [1E-2, 0.2E-2, 0.2E-2, 1, 32, 33, 8.3801E-6, 1.2662E-09, 0.05E-2, 4E+5, 30]; 

Hepatic = [6.6E-2, 0.22E-2,0.22E-2, 0, 0.36300E+10, 0.51251E-10, 9.3668E-6, 1.0505E-09, 0.049E-2, 4E+5, 31];

Gastric = [7.1E-2, 1.8E-2,1.8E-2, 0, 0.54100E+10, 0.34389E-10, 9.63070E-6, 8.7312E-10, 0.045E-2, 4E+5, 32];

Splenic = [6.3E-2, 0.275E-2,0.275E-2, 0, 0.23200E+10, 0.8019E-10, 8.8787E-6, 1.2487E-09, 0.054E-2, 4E+5, 33];

SuperiorMesenteric = [5.9E-2, 0.435E-2,0.435E-2, 0, 0.93000E+10, 0.20005E-09, 6.9676E-6, 1.7584-09, 0.069E-2, 4E+5, 34];

AbdominalAortaB = [1E-2, 0.6E-2,0.6E-2, 1, 36, 37, 3.095E-6, 2.5017E-09, 0.083E-2, 4E+5, 35]; 

LRenal = [3.2E-2, 0.260E-2,0.260E-2, 0, 0.11300E+10, 0.16464E-09, 8.9939E-6, 1.2077E-09, 0.052E-2, 4E+5, 36];

AbdominalAortaC = [1E-2, 0.59E-2,0.59E-2, 1, 38, 39, 3.5964E-6, 2.4148E-09, 0.083E-2, 4E+5,37];

RRenal = [3.2E-2, 0.26E-2,0.26E-2, 0, 0.11300E+10, 0.16464E-09, 8.9939E-6, 1.2077E-09, 0.052E-2, 4E+5, 38];

AbdominalAortaD = [10.6E-2, 0.58E-2, 0.548E-2, 1, 40, 41, 5.5958E-6, 2.045E-09, 0.082E-2, 4E+5, 39];

InferirorMesentric = [5.0E-2, 0.16E-2,0.16E-2, 0, 0.68800E+10, 0.27041E-10, 9.6873E-6, 7.7582E-10, 0.043E-2, 4E+5, 40];

AbdominalAortaeE = [1E-2, 0.52E-2,0.52E-2, 1, 42, 43, -0.2592E-5, 3.3951E-09, 0.078E-2, 4E+5, 41];

RCommonIlliac = [5.8E-2, 0.368E-2,0.35E-2, 1, 44, 45, 9.68951E-6, 7.59771E-10, 0.063E-2, 4E+5, 42];

LCommonIlliac = [5.8E-2, 0.368E-2,0.35E-2, 1, 50, 51, 9.68951E-6, 7.59771E-10, 0.063E-2, 4E+5, 43];

LExternalIliac = [14.4E-2, 0.32E-2,0.27E-2, 1, 46, 47, -0.64348E-6, 3.10359E-09, 0.055E-2, 4E+5, 44];

LInternalIliac = [5E-2, 0.2E-2,0.2E-2, 0, 0.79360E+10, 0.23443E-10, -0.18647E-4, 5.51695E-09, 0.050E-2, 4E+5, 45]; 

LFemoral = [44.3E-2, 0.259E-2,0.19E-2, 1, 48, 49, 9.68610E-6, 7.2179E-10, 0.052E-2, 4E+5, 46];

LDeepFemoral = [12.6E-2, 0.255E-2,0.186E-2, 0, 0.47700E+10, 0.39003E-10, 4.8860E-6, 6.569E-11, 0.046E-2, 16E+5, 47];

LPosteriorTibial = [32.1E-2, 0.247E-2,0.141E-2, 0, 0.47700E+10, 0.39003E-10, 4.6538E-6, 5.8506E-11, 0.051E-2, 16E+5, 48];

LAnteriorTibial = [34.3E-2, 0.13E-2,0.13E-2, 0, 0.55900E+10, 0.33281E-10, 4.072E-6, 4.2755E-11, 0.039E-2, 16E+5, 49];

RExternalIliac = [14.4E-2, 0.32E-2,0.27E-2, 1, 52, 53, -0.64348E-6, 3.10359E-09, 0.055E-2, 4E+5, 50];

RInternalIliac = [5E-2, 0.2E-2,0.2E-2, 0, 0.79360E+10, 0.23443E-10, -0.18647E-4, 5.51695E-09, 0.055E-2, 4E+5, 51]; 

RFemoral = [44.3E-2, 0.259E-2,0.19E-2, 1, 54, 55, 9.6861E-6, 7.2179E-10, 0.052E-2, 4E+5, 52];

RDeepFemoral = [12.6E-2, 0.255E-2,0.186E-2, 0, 0.47700E+10, 0.39003E-10, 4.88360E-6, 6.569E-11, 0.046E-2, 16E+5, 53];

RPosteriorTibial = [32.1E-2, 0.247E-2,0.141E-2, 0, 0.47700E+10, 0.39003E-10, 4.6538E-6, 5.8506E-11, 0.051E-2, 16E+5, 54];

RAnteriorTibial = [34.3E-2, 0.13E-2, 0.13E-2, 0, 0.55900E+10, 0.33281E-10, 4.072E-6, 4.2755E-11, 0.039E-2, 16E+5, 55];






System = [
    AscendingAorta;
    AorticArchA;
    Innominate;
    RSubclavianA;
    RCartoid;
    RVertebral;
    RSubclavianB;
    RRadial;
    RUlnarA;
    RInterosseus;
    RUlnarB;
    RInternalCartoid;
    RExternalCartoid;
    AorticArchB;
    LCartoid;
    LInternalCartoid;
    LExternalCartoid;
    ThoraicAortaA;
    LSubclavianA;
    Vertebral;
    LSubclavianB;
    LRadial;
    LUlnarA;
    LInterosseous;
    LUlnarB;
    Intercostals;
    ThoraicAortaB;
    AbdominalAortaA;
    CeliacA;
    CeliacB;
    Hepatic;
    Gastric;
    Splenic;
    SuperiorMesenteric;
    AbdominalAortaB;
    LRenal;
    AbdominalAortaC;
    RRenal;
    AbdominalAortaD;
    InferirorMesentric;
    AbdominalAortaeE;
    RCommonIlliac;
    LCommonIlliac;
    LExternalIliac;
    LInternalIliac;
    LFemoral;
    LDeepFemoral;
    LPosteriorTibial;
    LAnteriorTibial;
    RExternalIliac;
    RInternalIliac;
    RFemoral ;
    RDeepFemoral ;
    RPosteriorTibial ;
    RAnteriorTibial;
    ];


% Labelling the terminal nodes so we can first apply the boundary
% conditions here to find our relation between the reflection coefficients
% and indcidence coefficients for each artery. We then step back from
% this point and use the fact mass flux is conserved and pressure is
% continuous with the relations to find the relation between the incidence
% and reflection coefficients of the parent nodes


TerminalNodes = [6, 8, 10, 11, 12, 13, 16, 17, 20, 22, 24, 25, 26, 31, 32, 33, 34, 36, 38, 40, 45, 47, 48, 49, 51, 53, 54, 55];
% Matrix where each entry gives a junction with the first index giving the
% parent node and the latter two giving the daughter nodes.
             
MergingJunctions = [
    1,2,3;
    2,14,15;
    3,4,5;
    4,6,7;
    5,12,13;
    7,8,9;
    9,10,11;
    14,18,19;
    15,16,17;
    18,26,27;
    19,20,21;
    21,22,23;
    23,24,25;
    27,28,29;
    28,34,35;
    29,30,31;
    30,32,33;
    35,36,37;
    37,38,39;
    39,40,41;
    41,42,43;
    42,44,45;
    43,50,51;
    44,46,47;
    46,48,49;
    50,52,53;
    52,54,55;
    ];
    
            

                 
% Finding how big the arterial system is 
[N, NumVar] = size(System);

% Initialising the data
k = 10;
rho = 0.105E+4;
c0 = zeros(N,1);

% Code to create the PWV
for i = 1: N
   Artery = System(i,:);
   r1 = Artery(2);
   r2 = Artery(3);
   r = (r1+r2)/2;
   E = Artery(10);
   h = Artery(9);
   c0(i) = sqrt(E*h/2/rho/r);
   c0(i) = 1;
end


% Matrix which shows the relation between the incidence and reflection
% coefficient for a particular artery. (i,j) denotes node i in artery j.
RelationalCoefficients = zeros(k+1, N);


% Finding the relation between the coefficients at the terminal branches
% using the fact 2 windkessel model function.




% This is the code that loops through the terminal nodes to work out the
% relation

% Currently a capictor and resistor windkessel model. If we wish to make it
% just a resistor then set R = RT instead of the function of C and RT


for i = 1 : length(TerminalNodes)
    node = TerminalNodes(i);
    Artery = System(node,:);
    L = Artery(1);
    r = (Artery(2)+Artery(3))/2;
    A = pi*r^2;
    
    %L = Artery(1);
    %r1 = Artery(2);
    %r2 = Artery(3);
    %A = pi/2*(r1^2+r2^2);
    c = c0(node);
    Y = A / rho / c;
    CT = Artery(6);
    RT = Artery(5);
    R1 = 0.2*RT;
    R2 = RT -R1;
    RelationalCoefficients(1,node) = (RT * Y - 1) / (RT*Y +1);
    for j = 2 : k+1
        Qcoeff = (1i * (j-1) + (1 + R1/R2)/R1/CT);
        Pcoeff = (1i*(j-1)/R1) + (1/R1/R2/CT);
        R = (1i * (j-1) + (1 +R1/R2)/R1/CT)/((1i*(j-1)/R1) + (1/R1/R2/CT));
        constant = (R * Y - 1) / (R*Y +1);
        RelationalCoefficients(j,node) = constant * exp(1i * (j-1)*-2* L/c);
    end 
end



% The number of junctions we are considering.
Junctions = size(MergingJunctions,1);


% Code to use the relation between the incidence and relation coefficients
% of daugther nodes to then work out the relation in the parental nodes.
% Start from the end and work our way backwards

for i = Junctions:-1:1
    Art1 = System(MergingJunctions(i,1),:);
    Art2 = System(MergingJunctions(i,2),:);
    Art3 = System(MergingJunctions(i,3),:);
    L1 = Art1(1);
    R1 = (Art1(2)+Art1(3))/2;
    A1 = pi*R1^2;
    c1 = c0(MergingJunctions(i,1));
    Y1 = A1 / rho / c1;
    
    R2 = (Art2(2)+Art2(3))/2;
    A2 = pi*R2^2;
    c2 = c0(MergingJunctions(i,2));
    Y2 = A2 / rho / c2;
    
    R3 = (Art3(2)+Art3(3))/2;
    A3 = pi*R3^2;
    c3 = c0(MergingJunctions(i,3));
    Y3 = A3 / rho / c3;
    
    
    
   
    for j = 1 : k+1
        alpha = RelationalCoefficients(j, MergingJunctions(i,2));
        beta = RelationalCoefficients(j, MergingJunctions(i,3));
        f2 = (1-alpha)/(1+alpha);
        f3 = (1-beta)/(1+beta);
        Quotient = Y1 + Y2*f2 + Y3 * f3;
        RelationalCoefficients(j,MergingJunctions(i,1)) =     exp(1i*-2*(j-1)*L1/c1) * (Y1 - Y2*f2 - Y3*f3)/Quotient ; 
    end
end


%Incidence and reflection matrices where (i,j) denotes node i in artery j.
Incidence = zeros(k+1,N);
Reflection = zeros(k+1,N);
Incidence(:,1) = ones(k+1,1);


% We then work forwards using our relational coefficients and the values of
% the parental nodes to calculate the incidence and reflection coefficients
% of the daughter nodes.

for i = 1 : Junctions
    Art1 = System(MergingJunctions(i,1),:);
    Art2 = System(MergingJunctions(i,2),:);
    Art3 = System(MergingJunctions(i,3),:);
    c = c0(MergingJunctions(i,1));
    L1 = Art1(1);
    for j = 1:k+1
        alpha = RelationalCoefficients(j, MergingJunctions(i,2));
        beta = RelationalCoefficients(j, MergingJunctions(i,3));
        Reflection(j,MergingJunctions(i,1)) = Incidence(j,MergingJunctions(i,1)) * RelationalCoefficients(j,MergingJunctions(i,1));
        Incidence(j,MergingJunctions(i,2)) = 1/(alpha+1) * (Incidence(j,MergingJunctions(i,1))* exp(-1i*(j-1)*L1/c) + Reflection(j,MergingJunctions(i,1))* exp(1i*(j-1)*L1/c));
        Incidence(j,MergingJunctions(i,3)) = 1/(beta+1) * (Incidence(j,MergingJunctions(i,1))* exp(-1i*(j-1)*L1/c) + Reflection(j,MergingJunctions(i,1))* exp(1i*(j-1)*L1/c));
    end
end


% Finally use the relations to find the values of them at the end

for i = 1 : length(TerminalNodes)
    node = TerminalNodes(i);
    for j = 1 : k+1
        Reflection(j,node) = RelationalCoefficients(j,node) * (Incidence(j, node));
    end 
end


% Code to turn real fourier constants for sin and cos into comple for e^ix.
% First index indicates whether sin or cos
% Second index highlights the fourier node
fourierCs = zeros(2,k+1);
fourierCs(1,1) = 0.86393E-4;
fourierCs(1,2) = -0.88455E-4;
fourierCs(1,3) = -0.52515E-4;
fourierCs(1,4) = 0.86471E-4;
fourierCs(1,5) = -0.26395E-4;
fourierCs(1,6) = -0.12987E-4;
fourierCs(1,7) = 0.20133E-5;
fourierCs(1,8) = 0.70896E-5;
fourierCs(1,9) = 0.32577E-5;
fourierCs(1,10) = -0.56573E-5;
fourierCs(1,11) = -0.19302E-5;

fourierCs(2,2) = 0.13368E-3;
fourierCs(2,3) = -0.12280E-3;
fourierCs(2,4) = 0.22459E-4;
fourierCs(2,5) = 0.22693E-4;
fourierCs(2,6) = 0.22398E-5;
fourierCs(2,7) = -0.22315E-4;
fourierCs(2,8) = 0.10065E-4;
fourierCs(2,9) = -0.21066E-5;
fourierCs(2,10) = 0.90633E-5;
fourierCs(2,11) = -0.85422E-5;


complexFvals = zeros(1, k+1);
complexFvals(1)= fourierCs(1,1);
for h = 2 : k+1
    alpha = fourierCs(1,h);
    beta = fourierCs(2,h);
    complexFvals(h) = alpha - 1i *beta ;
end

% Finding what the values of the fourier constants have to be so that P(t)
% = f(t), our intial flux field at the moment.

scaledcomplexFvals = zeros(1,k+1);





for h = 1 : k+1
    Art1 = System(1,:);
    L1 = Art1(1);
    r1 = Art1(2);
    r2 = Art1(3);
    r = (r1+r2)/2;
    A1 = pi*r^2;
    c=c0(1);
    Y1 = A1 / rho /c;
    scaledcomplexFvals(h) = complexFvals(h) * ( 1 /(Y1*(Incidence(h,1)  - Reflection(h,1))));
end



% Creating the mesh of t vals and xvals.
xsteps = 1000;
tvals = 0:0.1:10;
tsteps = length(tvals);

%Creating the tensors to hold the information on flux and pressure.
Q = zeros(N, tsteps,xsteps);
P = zeros(N, tsteps,xsteps);
forward = zeros(N,tsteps,xsteps);
backward = zeros(N,tsteps,xsteps);


% Code that goes through each artery and calculates the flux and pressure
% for each.


for i = 1: N
    Artery = System(i,:);
    L = Artery(1);
    r1 = Artery(2);
    r2 = Artery(3);
    r=(r1+r2)/2;
    A = pi*r^2;
    c=c0(i);
    Y = A / rho / c;   
    xvals = linspace(0, L, xsteps);
    for j = 1: tsteps
        for h = 1: xsteps
            tn = tvals(j);
            xm = xvals(h);
            p=0;
            q=0;
            forward1=0;
            backward1=0;
            f = fourierVec(k,tn - xm/c,scaledcomplexFvals);
            b = fourierVec(k,tn +xm/c,scaledcomplexFvals);
            
            for m = 1 : length(f)
                q = q + Y *(Incidence(m,i)*f(m) - Reflection(m,i)*b(m));
                p = p + (Incidence(m,i)*f(m) + Reflection(m,i)*b(m));
                
            end
            for m = 1 : length(f)
                forward1 = forward1 + (Incidence(m,i)*f(m));
                backward1 = backward1 + (Reflection(m,i)*b(m));    
            end
            
            forward(i,j,h) = Y *real(forward1);
            backward(i,j,h) = Y *real(backward1);
            
            Q(i,j,h) = real(q);
            P(i,j,h) = real(p);
        end
    end
end
%%

%Plotting pressure and flux over time at various points of the arterial
%system


figure(1)
plot(tvals/2/pi*1.1,P(1,:,1))
hold on 
plot(tvals/2/pi*1.1,P(52,:,1))
plot(tvals/2/pi*1.1,P(28,:,1))
plot(tvals/2/pi*1.1,P(22,:,1))
hold off
legend('Ascending Aorta','Femoral Artery','Abdominal Aorta','Radial')
xlabel('time')
ylabel('Pressure, Pa ')
title('Pressure during over time at various points in the arterial system')

figure(2)
plot(tvals/2/pi*1.1,1000000*Q(1,:,1))
hold on 
plot(tvals/2/pi*1.1,1000000*Q(52,:,1))
plot(tvals/2/pi*1.1,1000000*Q(28,:,1))
plot(tvals/2/pi*1.1,1000000*Q(22,:,1))
hold off
legend('Ascending Aorta','Femoral Artery','Abdominal Aorta','Radial')
xlabel('time')
ylabel('Flux ml/s')
title('Flux during over time at various points in the arterial system')

%%


%Comparing the PWV to the lengths in the system
min(c0)
max(c0)
min(System(:,1))
max(System(:,1))

%%
for i = 1:Junctions
    figure(i)
    plot(tvals,Q(MergingJunctions(i,1)))
    hold on
    plot(tvals,Q(MergingJunctions(i,2))+Q(MergingJunctions(i,3)))
    hold off
end





%%

% Test to ensure that pressure is continuous and flux is conserved at
% junctions.
for i = 1:Junctions
    i;
    max(P(MergingJunctions(i,1),:,end)-P(MergingJunctions(i,2),:,1))
    max(P(MergingJunctions(i,1),:,end)-P(MergingJunctions(i,3),:,1))
    m = max(Q(MergingJunctions(i,1),:,end)-(Q(MergingJunctions(i,2),:,1)+Q(MergingJunctions(i,3),:,1)));
end

%%

y2=zeros(tsteps,1);

for i = 1: tsteps
    y2(i) = real(sum(fourierVec(k,tvals(i),complexFvals)));
end

% Shows that our flux field matches our inital in flow field.
figure(2)
plot(tvals'/2/pi*1.1, 1000000*y2);
xlabel('time (s)')
ylabel('Flux, ml/s ')
title('Flow at Ascending Aorta inlet')



%%
figure(1)
plot(tvals,forward(1,:,end))
hold on
plot(tvals,backward(1,:,end))
plot(tvals,forward(2,:,1))
plot(tvals, backward(2,:,1))
plot(tvals,forward(3,:,1))
plot(tvals, backward(3,:,1))
hold off
title('Showing the forward and backward running waves at the first merging junction')
legend('Forward Wave Ascending Aorta', 'Backward Wave Ascending Aorta','Forward Wave Aortic Arch A','Backward Wave Aortic Arch A' ,'Forward Wave Innominate','Backward Wave Innominate')


