L1 = 42;
L2=32;
L3=52.8;

F1=20000;
F2=25000;

MatsSets(1).E = 5.2e6;
MatsSets(1).A = 0.8;
%MatsSets(1).rho = 7.1;

%MatsSets(2).E = 70e9;
%MatsSets(2).A = 3e-3;
%MatsSets(2).rho = 3.5;

PD.N = 4;
PD.NodePos = [0, 0, 0;
              L1, 0, 0;
              L1, L2, 0;
              0, L2, 0];
PD.NE = 4;
PD.ElmConnect = [1, 2;
                 2, 3;
                 3, 4;
                 3, 1]
PD.NM = 1;
PD.MatsSets = MatsSets;
PD.ElmMats = [1;
              1;
              1;
              1]
PD.BCType = [1, 1, 1;
             0, 0, 1;
             0, 0, 1;
             1, 1, 1]
PD.BCVal = [0, 0, 0;
            F1, 0, 0;
            0, -F2,0;
            0, 0, 0];

PDans = PD_truss_static(PD);