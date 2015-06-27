clear;
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arial Cyr'); 
% k = 1.01:0.5:400;
% z0 = 0:0.5:50;

k = 0:0.005:1.5- 0.001;
z0 = 0:0.005:4*pi;

% k = 0:0.01:1.5 - 0.001;
% z0 = 0:0.01:250;

% k = 0:0.01:2;
% z0 = 0:0.01:4*pi;

% k = 0:0.00005:0.01;
% z0 = 0:0.00005:0.1;

%l = 2 * 7.3 * 10^(-5);
%l = 0.5;

maxZ = k;
maxS = k;


minZ = k;
minS = k;

% testMaxZ = 1.9633*sqrt(k);
% testMaxK = k(testMaxZ(:) <= max(z0));
% testMaxZ = testMaxZ(testMaxZ(:) <= max(z0));
% testMaxS = atan((0.5135 + testMaxK.^0.5 + 0.4856 * testMaxK)./( 0.4856 + testMaxK.^0.5 + 0.5135 * testMaxK));
% 
% 
% 
% testMinZ = (3 * k.^2).^(1/3);
% testMinS = atan((9 * k).^(1/3) / 2);
% 
% 
 lokMaxZ = 3.5343*sqrt(k);
% lokMaxK = k(lokMaxZ(:) <= max(z0));
% lokMaxZ = lokMaxZ(lokMaxZ(:) <= max(z0));
% lokMaxS = atan( (0.4993 + lokMaxK.^0.5 + 0.5006 * lokMaxK)./( 0.5006 + lokMaxK.^0.5 + 0.4993 * lokMaxK) );


% testMaxZ = 1.9633*sqrt(k);
% testMaxK = k(testMaxZ(:) <= max(z0));
% testMaxZ = testMaxZ(testMaxZ(:) <= max(z0));
% testMaxS = atan((0.5135 + testMaxK.^0.5 + 0.4856 * testMaxK)./( 0.4856 + testMaxK.^0.5 + 0.5135 * testMaxK));


w = -(2.*tan(z0).*tanh(z0) - 2)./(tan(z0) - tanh(z0));

temp1 = 1./(0.5*w + 0.5*((w + 2.0).*(w - 2.0)).^(1/2)).^2;
temp2 = 1./(0.5*w - 0.5*((w + 2.0).*(w - 2.0)).^(1/2)).^2;




[Z0, K] = meshgrid(z0,k);

res = zeros(length(k),length(z0));
res2 = zeros(length(k),length(z0));
c1 = zeros(length(k),length(z0));
c2 = zeros(length(k),length(z0));

% for i = 1:length(z0)
%     for t = 1:length(k)
%         C = funcArc(z0(1,i), k(1,t));
%         %if C(3,1) + C(4,1) ~= 0
%             res(t,i) = atan((C(3,1) - C(4,1))/(C(3,1) + C(4,1))); 
%             c1(t,i) = C(3,1);
%             c2(t,i) = C(4,1);
%         %end
%     end
% end

for t = 1:length(k)
    max = 0;
    min = 10;
  
    for i = 1:length(z0)
%      a = 2*k(1,t) * cosh(z0(1,i)/k(1,t)^0.5) - (k(1,t)^0.5 - k(1,t)^1.5) * sin(z0(1,i)/k(1,t)^0.5) + (k(1,t)^0.5 + k(1,t)^1.5) * sinh(z0(1,i)/k(1,t)^0.5);
%      b = 2*k(1,t) * cosh(z0(1,i)/k(1,t)^0.5) + (k(1,t)^0.5 - k(1,t)^1.5) * sin(z0(1,i)/k(1,t)^0.5) + (k(1,t)^0.5 + k(1,t)^1.5) * sinh(z0(1,i)/k(1,t)^0.5);
        
%        p = z0(1,i) * ( l /  ( 0.5 * k(1,t) )  )^0.5;
%        
%        a = 2*k(1,t) * cosh(p) - (k(1,t)^0.5 - k(1,t)^1.5) * sin(p) + (k(1,t)^0.5 + k(1,t)^1.5) * sinh(p);
%        b = 2*k(1,t) * cosh(p) + (k(1,t)^0.5 - k(1,t)^1.5) * sin(p) + (k(1,t)^0.5 + k(1,t)^1.5) * sinh(p);
       
       
%        p = 2*z0(1,i)*k(1,t)^-0.5;
%        
%        a = 2*k(1,t)^0.5 * cosh(p) - (1 - k(1,t)) * sin(p) + (1 + k(1,t)) * sinh(p);
%        b = 2*k(1,t)^0.5 * cosh(p) + (1 - k(1,t)) * sin(p) + (1 + k(1,t)) * sinh(p);
        
        p = z0(1,i);
        a = 2*k(1,t)^0.5 * cosh(p) - (1 - k(1,t)) * sin(p) + (1 + k(1,t)) * sinh(p);
        b = 2*k(1,t)^0.5 * cosh(p) + (1 - k(1,t)) * sin(p) + (1 + k(1,t)) * sinh(p);

            
        res(t,i) = atan(a/b);
       
        if res(t,i) >= max
            max = res(t,i);
            maxZ(1,t) = z0(1,i);
            maxS(1,t) = max;
        end
        
        if res(t,i) <= min
            min = res(t,i);
            minZ(1,t) = z0(1,i);
            minS(1,t) = min;
        end
        
    end
end

figure;  hold on;
[c,h]=contour (Z0,K,res,0:0.02:0.784);  hold on;
[c,h]=contour (Z0,K,res,0.01:0.02:0.784);
[c,h]=contour (Z0,K,res,0.784:0.0001:0.79); 
[c,h]=contour (Z0,K,res,0.79:0.004:1.5);
[c,h]=contour (Z0,K,res,pi/4,'g','LineWidth',2); hold on;

% plot(minZ,k,'ro'); hold on;
% plot(maxZ,k,'bo'); hold on;

testLOL = zeros(4,length(k));
for i = 1:4
    % testLOL(i,:) = 0.5 * (pi/4 + pi * (i-1)) .* sqrt(k);
    testLOL(i,:) = 0.7884 + pi * (i-1); 
    % testLOL =  testLOL(testLOL(i,:) <= max(z0));
    plot(testLOL(i,:),k,'yo'); hold on;
end

% nevMax = abs(maxZ - testLOL(2,:));
% nevMin = abs(minZ - testLOL(1,:));
% nevLoc = abs(lokMaxZ - testLOL(3,:));
% axis([ 0 4 0 1- 0.001]);




% plot3(Z0,K,res); hold on;
% plot3(maxZ,k,maxS,'bo'); hold on;
% plot3(minZ,k,minS,'ro'); hold on;
% ylabel( 'k-/k+' ); xlabel( 'z0' ); zlabel( 'угол' );
% grid;


% figure;
% mesh (Z0,K,res); hold on;
% plot3(lokMaxZ,lokMaxK,lokMaxS,'yo');
% plot3(testMaxZ,testMaxK,testMaxS,'ro');
% plot3(maxZ,k,maxS,'bo'); hold on;

% plot3(minZ,k,minS,'ro'); hold on;
% plot3(minZ,k,minS,'ro'); hold on;







% figure;  hold on;
% [c,h]=contour (Z0,K,res,0:0.02:0.784); clabel(c,h); hold on;
% [c,h]=contour (Z0,K,res,0.01:0.02:0.784);
% [c,h]=contour (Z0,K,res,0.7851:0.001:0.79); clabel(c,h);
% [c,h]=contour (Z0,K,res,0.79:0.004:1.5);
% [c,h]=contour (Z0,K,res,pi/4,'g','LineWidth',2);
% 
% plot(minZ,k,'ro'); hold on;
% plot(maxZ,k,'go'); hold on;


% plot(maxZ,k,'bo'); hold on;
% % % plot(testMinZ,k,'go'); hold on;
% % % plot(testMinZ,k,'go');
% plot(minZ,k,'ro'); 
% plot(lokMaxZ,lokMaxK,'yo');
% 
% plot([0 5],[1 1],'g','LineWidth',3);
% plot([0 0],[0 1.5],'g','LineWidth',3);

%[c,h]=contour (Z0,K,res,0.798);clabel(c,h); hold on; [c,h]=contour (Z0,K,res,0.795); ;clabel(c,h);

% xlabel('z0','FontName','Arial Cyr','FontSize',25);
% ylabel('k-/k+','FontName','Arial Cyr','FontSize',25);
% zlabel('”гол','FontName','Arial Cyr','FontSize',25);

xlabel('$\eta$','FontName','Arial Cyr','FontSize',25,'Interpreter', 'latex');
ylabel('$\kappa$','FontName','Arial Cyr','FontSize',25,'Interpreter', 'latex');
zlabel('”гол','FontName','Arial Cyr','FontSize',25);


% x = log(k(1,2:length(k)));  'Interpreter', 'latex'
% y = log(minZ(1,2:length(k)));
% p = polyfit(x, y, 1);
% q2 = exp(p(1,2))*k.^(p(1,1));