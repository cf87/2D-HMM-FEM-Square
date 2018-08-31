function b=neumannbc(g, q)
%% Assemble boundary
 
if nargin<2, q='0'; end

b1      = nbc(q,g) ;
b2      = nbc(q,g); 
b3      = nbc(q,g); 
b4      = nbc(q,g); 



b       = makeBoundaryMatrix( b1 , b2 , b3 , b4 ,0);



end


function b=makeBoundaryMatrix(b1,b2,b3,b4,nd)

q=     strvcat(    b1.q,  b2.q,  b3.q,  b4.q);
g=     strvcat(    b1.g,  b2.g,  b3.g,  b4.g);

lq=length(q(1,:));
lg=length(g(1,:));

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


end

function bdry= nbc(q,g)

bdry.q=q;
bdry.g=g;
bdry.h= [];
bdry.r=[];

end

