function PD=hexagon_torsion_example(a)

PD.VertexList = [0, a;
                 a, 0;
                 a, -a;
                 0, -2*a;
                 -a, -a;
                 -a, 0;
                 0, a];
PD.InitEdgeLen = 0.65;
PD.BBox = [-3*a, -3*a;2*a, 3*a];
PD.RHS = -1.0;

PD=PD_torsion_poly(PD,1);
