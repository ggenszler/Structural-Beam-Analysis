% function solves the reaction forces for a simply supported beam
function rxns=statics(supports,load)

    [m,n]=size(load);
    pfy=0;
    dfy=0;
    pm=0;
    dm=0;
    
    for i=1:m
        if load(i,1)==load(i,3) && load(i,4)==1
            pfy=pfy-load(i,2);
            pm=pm-load(i,2)*(load(i,1)-supports(1,2));
            % calculates the amount of force due to a distributed load and
            % the moment caused by it
        elseif load(i,1)==load(i,3) && load(i,4)==0
            pfy=pfy;
            pm=pm+load(i,2);
            % accounts for an applied moment in the moment sum
        else
            dfy=dfy-load(i,2)*(load(i,3)-load(i,1));
            dm=dm-load(i,2)*(load(i,3)-load(i,1))*(1/2*(load(i,1)+load(i,3))-supports(1,2));
            % calculates the amount of force due to a point load and the
            % moment caused by it
        end
    end
    
    fy=pfy+dfy;
    mm=pm+dm;
    % combines all moments and loads
    
    if supports(1,1)==0 || supports(2,1)==0
        rxns=[fy, mm];
        % assigns the sums as force in the y direction and reactional
        % moment
    else
        By=(mm)/(supports(2,2)-supports(1,2));
        Ay=(fy-By);
        rxns=[Ay, By];
        % performs sum of moments about the left most support and uses the
        % answer to solve sum of forces in the y direction
    end
end