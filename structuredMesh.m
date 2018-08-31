function [p e t midpoints]= structuredMesh(x,y)
% Creates a simple structured rectanguar triangulated mesh of x \times y 
% according to matlabs inititmesh specifications


%number of triangles along a horizontal edge and a vertical edge
numx=length(x)-1;
numy=length(y)-1;

% In the Point matrix p, the first and second rows contain x- and
% y-coordinates of the points in the mesh.
[X,Y] = meshgrid(x,y);
X1=reshape(X',length(x)*length(y),1);
Y1=reshape(Y',length(x)*length(y),1);
p=[X1'; Y1']; 


% In the Triangle matrix t, the first three rows contain indices to the
% corner points, given in counter clockwise order, and the fourth row
% contains the subdomain number

t       = ones(4, 2*numx*numy);
midpoints = zeros(2,length(t));
ind     = 1:numx;
for j=1:numy
    v                   = (j-1)*(numx+1)+ (1:numx);
    t(1:3, ind)         = [v; v+numx+2; v+numx+1]; %top left triangles
    midpoints(1,ind)    = mean(reshape(p(1,t(1:3,ind)), [3, length(ind)]));
    midpoints(2,ind)    = mean(reshape(p(2,t(1:3,ind)), [3, length(ind)]));
    
    ind=ind+numx;
    t(1:3, ind)         = [v; v+1; v+numx+2]; %bottom right triangles
    midpoints(1,ind)    = mean(reshape(p(1,t(1:3,ind)), [3, length(ind)]));
    midpoints(2,ind)    = mean(reshape(p(2,t(1:3,ind)), [3, length(ind)]));
    ind                 = ind+numx;
end



% In the Edge matrix e, the first and second rows contain indices of
%the starting and ending point, the third and fourth rows contain the
%starting and ending parameter values, the fifth row contains the edge
%segment number, and the sixth and seventh row contain the left- and
%right-hand side subdomain numbers.
%

% fills in the edge segments with a vector v of starting points, v+inc
% ending points, and seg_num segment number
edge=@(v, inc, seg_num) [v; v+inc; zeros(size(v)); ones(size(v)); seg_num*ones(size(v)); zeros(size(v)); ones(size(v))];

s       = 1:numx;
south   = edge(s,1,1);
       
ea      = numx+1:(numx+1): (numy)*(numx+1);
east    = edge(ea, numx+1, 2);

n       = (numx+1)*(numy+1):-1: (numx+1)*numy+2;
north   = edge(n, -1, 3);

w       = ((numx+1)*numy+1):-1*(numx+1):numx+2;
west    = edge(w,-1*(numx+1),4);

e       = [south east north west];
end


