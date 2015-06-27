function [ asimpt1, asimpt2, asimpt3, asimpts, exact, sj7 ] = asimptotics(k, s, sorz)

 if sorz == 0
    asimpt1 = 0.5 * k .* ((k + k.^-1)/3).^(-1/3);
    asimpt2 = 0.5 * k .* (((k + k.^-1)/3).^(-1/3) - (k + k.^-1).^-1);
    asimpt3 = 0.5 * k .* (((k + k.^-1)/3).^(-1/3) - (k + k.^-1).^-1 + (18/(35*3^(1/3)))*(k + k.^-1).^-(5/3));
 else
    asimpt1 = ((k + k.^-1)/3).^(-1/3);
    asimpt2 = ((k + k.^-1)/3).^(-1/3) - (k + k.^-1).^-1;
    asimpt3 = ((k + k.^-1)/3).^(-1/3) - (k + k.^-1).^-1 + (18/(35*3^(1/3)))*(k + k.^-1).^-(5/3);
 end

fa = tan(s) - tanh(s);
sj = zeros(1,2);
counter = 1;
for i = 2:length(s)
    if fa(i)*fa(i-1) <= 0 && fa(i-1) <=0 && s(i) > 0.1
        sj(counter) = s(i);
        counter = counter + 1;
    end
end

sj1 = cos(sj).^2 ./ (cos(sj).^2 - cosh(sj).^2);

asimpts = zeros(length(sj),length(k));
temp = zeros(length(sj),length(k));

for t = 1:length(sj)
    for i = 1:length(k)
        if sorz == 0
            asimpts(t,i) = 0.5 * sj(t) * k(i) + 0.5 * sj1(t) * k(i)^3;
            temp(t,i) = 0.5 * sj(t) * k(i) ;
        else
            asimpts(t,i) = sj(t) + sj1(t)* k(i)^2;
            temp(t,i) = sj(t);
        end
    end
end

sj7 = temp;


exact = zeros(length(k),3);
f = zeros(length(k),length(s));

for i = 1:length(k)
    w = k(i) + k(i)^-1;
    f(i,:) = tan(s) - (2 + w * tanh(s))./(w + 2 * tanh(s));
    counter = 1;
    for t = 2:length(s)
        if f(i,t)*f(i,t-1) <= 0 && f(i,t-1) <=0
            exact(i,counter) = s(t);
            counter = counter + 1;
        end
    end
    %clear f
end
if sorz == 0
    for t = 1:length(exact(1,:))
        for i = 1:length(k)
            exact(i,t) = 0.5 * exact(i,t) * k(i);
        end
    end 
end


end

