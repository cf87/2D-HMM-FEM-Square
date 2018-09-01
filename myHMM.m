function [u, A_HMM] = myHMM(macro_mesh, micro_mesh, aeps,pde, opts)
% Disclaimer: this code was written in 2014 and it appears assema and
% assempde are now out of date.

% The microscopic cell problems satisfy $\langle nabla \phi^\epsilon
% \rangle_{I_\delta} = G$ for some fixed vector $G$

% The cell problems either have Dirichlet or Neumann BC
if ~isfield(opts, 'bc'), opts.bc='dirichlet'; end;
switch opts.bc
    case 'dirichlet'
        bc1     = dirichletbc(sprintf('%d*x+%d*y', 1, 0));
        bc2     = dirichletbc(sprintf('%d*x+%d*y', 0, 1));
    case 'neumann'
        
        bc1     = neumannbc(sprintf('%d*nx+%d*ny', 1, 0));
        bc2     = neumannbc(sprintf('%d*nx+%d*ny', 0, 1));
end

aeps_mp = aeps(micro_mesh.mp(1,:), micro_mesh.mp(2,:)); %multiscale coefficient a

% Cell problems (there is one in each dimension)

     [K,M,F]         = assema(micro_mesh.p,micro_mesh.t,aeps_mp,'0','0');
     [Q1,G1,H1,R1]   = assemb(bc1,micro_mesh.p,micro_mesh.e);
     [Q2,G2,H2,R2]   = assemb(bc2,micro_mesh.p,micro_mesh.e);
     
     phi1            = assempde(K,M,F,Q1,G1,H1,R1) ;
     phi2            = assempde(K,M,F,Q2,G2,H2,R2) ;

[ phi1x,phi1y ] = pdegrad(micro_mesh.p,micro_mesh.t,phi1);
[ phi2x,phi2y ] = pdegrad(micro_mesh.p,micro_mesh.t,phi2);
rs=@(x) reshape(x, [length(x)/length(macro_mesh.mp) length(macro_mesh.mp) ]);

[ agx_phi1,agy_phi1 ] = pdecgrad(micro_mesh.p,micro_mesh.t,aeps_mp,phi1);
[ agx_phi2,agy_phi2 ] = pdecgrad(micro_mesh.p,micro_mesh.t,aeps_mp,phi2);



% Given these cell solutions we want to solve 
% \langle a^\epsilon \nabla phi^\epsilon \rangle_{I_delta} = A_HMM \langle phi^\epsilon \rangle_{I_delta}%

% First, find cell averages
cell_avg    = @(x) mean(rs(x),1);
V1    = cell_avg(phi1x); V2=cell_avg(phi2x); V3=cell_avg(phi1y); V4=cell_avg(phi2y);
AV1   = cell_avg(agx_phi1); AV2=cell_avg(agx_phi2); AV3=cell_avg(agy_phi1); AV4=cell_avg(agy_phi2);

%Assemble them into a large linear system and solve for A_HMM
V       = blckdiag(V1, V2, V3, V4);
AV      = blckdiag(AV1, AV2, AV3, AV4);
A_HMM   = unblckdiag(AV/V);



%% Macroscopic solution
[K,M,F] = assema(macro_mesh.p,macro_mesh.t,A_HMM',0,pde.f);
u       = assempde(K,M,F,pde.Q,pde.G,pde.H,pde.R) ;

end

function VV=blckdiag(V1, V2, V3, V4)
row=@(x,N)x(1:N);
Vdiag=row([V1; V4], 2*length(V1));
Vp=row([V2; 0*V2], length(Vdiag)-1);
Vm=row([V3; 0*V3], length(Vdiag)-1);
VV=diag(Vdiag)+diag(Vp, 1)+diag(Vm, -1);
end

function V=unblckdiag(VV)
N=length(VV)/2;

Vdiag=diag(VV);

Vp=diag(VV,1);
Vm=diag(VV,-1);

V1=Vdiag(1:2:2*N);
V2=Vp(1:2:2*N);
V3=Vm(1:2:2*N); 
V4=Vdiag(2:2:2*N);



V=[V1 V2 V3 V4];
end
