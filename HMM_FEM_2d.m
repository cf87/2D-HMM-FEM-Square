

%% Quick and dirty HMM-FEM 2D solver
% Want to solve $\nabla \cdot (A_HMM \nabla U_HMM) = f$ on $\Omega = [0, 1]\times
% [0, 1]$, $U_HMM = g$ on $\partial \Omega$ (Dirichlet) or $\nabla U_HMM
% \cdot n = g$ (Neumann)

%% macroscopic RHS d and BC rhs g
pde.g           = '0'; 
pde.f           = '1';


%% a^epsilon and HMM options
epsilon     = 1/100;
a           = @(x1, x2, y1, y2) 1.1 + sin(2*pi*y2);
aeps        = @(x,y) a(x,y, x/epsilon, y/epsilon);



cellSize    =  2;
nMicro      = 30;
nMacro                  = 15;

load('g'); % I used pdetool to generate a square and saved it as the file 'g.mat';
hmmopts     = struct('epsilon', epsilon, 'nMicro',nMicro, 'deltax',cellSize*epsilon,'deltay',cellSize*epsilon, 'cellSize', cellSize, bc, 'Dirichlet');
%% This takes forever, so do it one time and save the mesh as a .mat file.
[macro_mesh, micro_mesh]=hmmmesh(nMacro,g, hmmopts);

figure;
pdemesh(micro_mesh.p, micro_mesh.e, micro_mesh.t); hold on;
pdemesh(macro_mesh.p, macro_mesh.e, macro_mesh.t);

%% Macroscopic boundary conditions
pde.bc='dirichlet';
switch pde.bc
    case 'dirichlet'
        b     = dirichletbc(pde.g);
    case 'neumann'
        b     = neumannbc(pde.g);
end

[pde.Q, pde.G, pde.H, pde.R] = assemb(b,macro_mesh.p, macro_mesh.e);




%% Produce HMM solution and coefficent A_HMM at the macro_mesh midpoints
[U_HMM, A_HMM] = myHMM(macro_mesh, micro_mesh, aeps,pde, hmmopts);
figure;
pdesurf(macro_mesh.p,macro_mesh.t,U_HMM);




