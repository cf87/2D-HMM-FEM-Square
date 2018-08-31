

%% Quick and dirty square macroscopic mesh


%%
pde.g                = '0';
pde.sigma          = @(x,y,l)0*x;
pde.K                 = 1;
pde.alpha             = 1;
pde.omega           = 0;
pde.alpha       =ones(1, pde.K);
pde.rhs = '1';

%% a^epsilon and HMM optionz
epsilon     = 1/100;
a           = @(x1, x2, y1, y2) 1.1 + sin(2*pi*y2);
aeps        = @(x,y) a(x,y, x/epsilon, y/epsilon);



cellSize    =  2;
nMicro      = 30;
nMacro                  = 15;

load('g'); % I used pdetool to generate a square and saved it as the file 'g.mat';
hmmopts     = struct('epsilon', epsilon, 'nMicro',nMicro, 'deltax',cellSize*epsilon,'deltay',cellSize*epsilon, 'cellSize', cellSize);
%%
[macro_mesh, micro_mesh]=hmmmesh(nMacro,g, hmmopts);
%%
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




%% Solve HMM
[u, A_HMM] = myHMM(macro_mesh, micro_mesh, aeps,pde, hmmopts);
figure;
pdesurf(macro_mesh.p,macro_mesh.t,u);




