% solves the displacement at a given distance along the beam
function d=displacement(y,rxns,load,supports,l,E,I,units)

[m1,n1]=size(supports);
[m2,n2]=size(load);

syms x;
dcoord=[];
tcoord=[];

for i=1:m1
    if supports(i,1)==0
        dcoord=[dcoord;
                supports(i,2)    rxns(1)*units/6*x^3-rxns(2)*10^3/2*x^2];
        tcoord=[tcoord;
                supports(i,2)    rxns(1)*units/2*x^2-rxns(2)*10^3*x];
    else
        dcoord=[dcoord;
                supports(i,2)    rxns(i)*units/6*(x-supports(i,2))^3];
    end
end

for i=1:m2
    if load(i,1)==load(i,3) && load(i,4)==1
        dcoord=[dcoord;
                load(i,1)    load(i,2)*units/6*(x-load(i,1))^3];
    elseif load(i,1)==load(i,3) && load(i,4)==0
        dcoord=[dcoord;
                load(i,1)    -load(i,2)*units/2*(x-load(i,1))^2];
    else
        dcoord=[dcoord;
                load(i,1)    load(i,2)*units/24*(x-load(i,1))^4;
                load(i,3)    -load(i,2)*units/24*(x-load(i,3))^4];
    end
end

if size(tcoord,1)==0
    ;
else
    for i=1:m2
        if load(i,1)==load(i,3) && load(i,4)==1
            tcoord=[tcoord;
                    load(i,1)    load(i,2)*units/2*(x-load(i,1))^2];
        elseif load(i,1)==load(i,3) && load(i,4)==0
            tcoord=[tcoord;
                    load(i,1)    -load(i,2)*units*(x-load(i,1))];
        else
            tcoord=[tcoord;
                    load(i,1)    load(i,2)*units/6*(x-load(i,1))^3;
                    load(i,3)    -load(i,2)*units/6*(x-load(i,3))^3];
        end
    end
    tc(x)=tcoord;
end
% creates discontinuity function component matrices. Each application or
% ending of a load is assigned a distance and discontinuity term per row.
% This is all done in terms of the variable x. The first for loop is for
% supports. The second loop is for loading on a simply supported beam. The
% third loop is for loading on a cantilever beam.

[m3,n3]=size(dcoord);
[m4,n4]=size(tcoord);
dc(x)=dcoord;
d=0;
t=0;
A=0;
B=0;

syms c1 c2

Amag=dc(supports(1,2));

for i=1:m3
    if supports(1,2)>=dcoord(i,1)
        A=A+Amag(i,2);
    else
        continue
    end
end
A=A+c1*supports(1,2)+c2;

if supports(1,1)==0
    tmag=tc(supports(1,2));
    for i=1:m4
        if supports(1,2)>=tcoord(i,1)
            B=B+tmag(i,2);
        else
            continue
        end
    end
    B=B+c1;
else
    Bmag=dc(supports(2,2));
    for i=1:m3
        if supports(2,2)>=dcoord(i,1)
            B=B+Bmag(i,2);
        else
            continue
        end
    end
    B=B+c1*supports(2,2)+c2; 
end
% The for loops above create the necessary discontinuity function system in
% order to solve for the two constants from double integration. 
   
eqs=[A==0, B==0];
vars=[c1 c2];
    
[sol1, sol2]=solve(eqs, vars);
% For a simply supported beam, the displacements at the locations of the 
% supports will be zero. For a cantilever beam, displacement and rotation
% at the fixed support will be zero.
  
for i=1:m3
    if y>=dcoord(i,1)
        dmag=dc(y);
        d=d+dmag(i,2);
    else
        continue
    end
end
% With the constants now solved for, the displacement at the location in
% question can now be calculated. 
d=1/(E*I)*(d+sol1*y+sol2);
end