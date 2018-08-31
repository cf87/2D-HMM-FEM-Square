function [macro_mesh, micro_mesh]=hmmmesh(n,g, opts)


    [p e t]=poimesh(g,n,n);
  
    x = (p(1,t(1,:))+p(1, t(2,:))+p(1,t(3,:)) )/3;
    y = (p(2,t(1,:))+p(2, t(2,:))+p(2,t(3,:)) )/3;
    
    mp = [x;y];
    macro_mesh.p=p;macro_mesh.e=e;macro_mesh.t=t;
    macro_mesh.mp=mp;
    
    %% This takes a long time, it is filling in local meshes
    for k=1:length(mp);
        
        cx  = mp(1,k);
        cy  = mp(2,k);
        x   = (cx-opts.deltax/2):opts.deltax/opts.nMicro:(cx+opts.deltax/2);
        y   = (cy-opts.deltax/2):opts.deltay/opts.nMicro:(cy+opts.deltay/2);
        
        [sgp, sge, sgt, sgmp]= structuredMesh(x,y);
        
        if k==1;
            meshp = sgp;
            meshe = sge;
            mesht = sgt;
            meshmp = sgmp;
            l=length(sgp);
        else
            
            
            sge(7,:)= sge(7,:)+4*(k-1);
            sgt(1,:)= sgt(1,:)+(k-1)*l;
            sgt(2,:)= sgt(2,:)+(k-1)*l;
            sgt(3,:)= sgt(3,:)+(k-1)*l;
            sge(7,:)=k;
            sgt(4,:)=k;
            sge(1,:) =  sge(1,:)+(k-1)*l;
            sge(2,:) =  sge(2,:)+(k-1)*l;
            
            meshp= [meshp  sgp] ;
            meshe= [meshe  sge] ;
            mesht= [mesht  sgt] ;
            meshmp= [meshmp  sgmp] ;
            
            
            clear sgp sge sgt sgmp k ;
        end
         micro_mesh.p=meshp; micro_mesh.e=meshe; micro_mesh.t=mesht; micro_mesh.mp=meshmp;
        
    end
    xx=0:opts.dx:1;
    x=meshgrid(xx);y=x';
   





end



