function b=dirichletbc(g)
%% Create boundary for [0 1]\times [0 1] boundary with four segments
%% Assemble boundary
%just based on my test example, need edge info to do this better%
b1      = dbc( '1', g     ); 
b2      = dbc( '1', g    ); 
b3      = dbc( '1', g     ); 
b4      = dbc( '1', g    ); 


b       = makeBoundaryMatrix( b1 , b2 , b3 , b4 ,1);



end


function b=makeBoundaryMatrix(b1,b2,b3,b4,nd)

q=     strvcat(    b1.q,  b2.q,  b3.q,  b4.q);
g=     strvcat(    b1.g,  b2.g,  b3.g,  b4.g);
h=     strvcat(    b1.h,  b2.h,  b3.h,  b4.h);
r=     strvcat(    b1.r,  b2.r,  b3.r,  b4.r);

lq=length(q(1,:));
lg=length(g(1,:));
lh=length(h(1,:));
lr=length(r(1,:));


b1      =[1     nd   lq     lg     lh     lr     q(1,:)     g(1,:)  h(1,:) r(1,:)]';
b2      =[1     nd   lq     lg     lh     lr     q(2,:)     g(2,:)  h(2,:) r(2,:)]';
b3      =[1     nd   lq     lg     lh     lr     q(3,:)     g(3,:)  h(3,:) r(3,:)]';
b4      =[1     nd   lq     lg     lh     lr     q(4,:)     g(4,:)  h(4,:) r(4,:)]';

if nd==0
    b1      =[1     nd   lq     lg         q(1,:)     g(1,:)  ]';
    b2      =[1     nd   lq     lg        q(2,:)     g(2,:)  ]';
    b3      =[1     nd   lq     lg        q(3,:)     g(3,:)  ]';
    b4      =[1     nd   lq     lg        q(4,:)     g(4,:)  ]';
    
end

b=[b1 b2 b3 b4];


end
function bdry= dbc(h,r)

bdry.q='0';
bdry.g='0';
bdry.h= h;
bdry.r=r;

end

function bdry= nbc(q,r)

bdry.q=q;
bdry.g=r;
bdry.h= '0';
bdry.r='0';

end

