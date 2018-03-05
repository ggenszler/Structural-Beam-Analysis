% function graphs shear and moment diagrams using discontinuity functions
function [v,m]=shearmoment(rxns,load,supports,l,lstep)

% rxns=[ -523/108, -199/36];
% load =[ 2,        6, 2, 0;
% 2,       -2, 3, 1;
%  3,       -8, 4, 1;
%  3, 1603/108, 3, 1];
% supports=[0 0];
% l=4;
% lstep=.1;




    [m1,n1]=size(supports);
    [m2,n2]=size(load);

    syms x;
    vcoord=[];
    mcoord=[];
    
    % supports: generates correct discontinuity term from supports
    for i=1:m1
        if supports(i,1)==0
            vcoord=[vcoord;
                    supports(i,2)    rxns(1)];
            mcoord=[mcoord;
                    supports(i,2)    rxns(1)*(x-supports(i,2))-rxns(2)];
        else
            vcoord=[vcoord;
                    supports(i,2)    rxns(i)];
            mcoord=[mcoord;
                    supports(i,2)    rxns(i)*(x-supports(i,2))];
        end
    end
    
    % loading: generates correct discontinuity term from loading
    for i=1:m2
        if load(i,1)==load(i,3)
            vcoord=[vcoord;
                    load(i,1)    load(i,2)+x*0];
            mcoord=[mcoord;
                    load(i,1)    load(i,2)*(x-load(i,1))];
        else
            vcoord=[vcoord;
                    load(i,1)    load(i,2)*(x-load(i,1));
                    load(i,3)    -load(i,2)*(x-load(i,3))];
            mcoord=[mcoord;
                    load(i,1)    1/2*load(i,2)*(x-load(i,1))^2;
                    load(i,3)    -1/2*load(i,2)*(x-load(i,3))^2];
        end
    end
    
    [m3,n3]=size(vcoord);
    [m4,n4]=size(mcoord);
    
    vc(x)=vcoord;
    mc(x)=mcoord;
    
    v=[];
    m=[];
    
    % computes shear and moment values and stores them in arrays
    lcount=0;
    for y=0:l/lstep
        vnew=0;
        for i=1:m3
            if lcount>=vcoord(i,1)
                vmag=vc(lcount);
                vnew=vnew+vmag(i,2);
            else
                continue
            end
        end
        v=[v vnew];
        lcount=lcount+lstep;
    end
        
    lcount=0;
    for y=0:l/lstep
        mnew=0;
        for i=1:m4
            if lcount>=mcoord(i,1)
                mmag=mc(lcount);
                mnew=mnew+mmag(i,2);
            else
                continue
            end
        end
        m=[m mnew];
        lcount=lcount+lstep;
    end
end