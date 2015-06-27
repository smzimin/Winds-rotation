function [ res, minZ, minS, maxZ, maxS ] = result(k, z0, sorz)

res = zeros(length(k),length(z0));

maxZ = k;
maxS = k;

minZ = k;
minS = k;

for t = 1:length(k)
    max = 0;
    min = 10; 
    for i = 1:length(z0)
        
        if sorz == 0
           p = 2*z0(1,i)*k(1,t)^-1;
           a = 2*k(1,t) * cosh(p) - (1 - k(1,t)^2) * sin(p) + (1 + k(1,t)^2) * sinh(p);
           b = 2*k(1,t) * cosh(p) + (1 - k(1,t)^2) * sin(p) + (1 + k(1,t)^2) * sinh(p);
        else
            p = z0(1,i);
            a = 2*k(1,t) * cosh(p) - (1 - k(1,t)^2) * sin(p) + (1 + k(1,t)^2) * sinh(p);
            b = 2*k(1,t) * cosh(p) + (1 - k(1,t)^2) * sin(p) + (1 + k(1,t)^2) * sinh(p);
        end
        res(t,i) = atan(a/b);
       
        if res(t,i) >= max
            max = res(t,i);
            maxZ(1,t) = z0(1,i);
            maxS(1,t) = max;
        end
        
        if res(t,i) <= min
            min = res(t,i);
            minZ(t) = z0(i);
            minS(t) = min;
        end
        
    end
end

end

