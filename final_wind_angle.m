clear;
fontsz = 16;

per  = 180 / pi;

for risunok = 6
close
set(0,'DefaultAxesFontSize',fontsz,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',fontsz,'DefaultTextFontName','Arial Cyr'); 
% выбор рисунка (в каких пределах и с какой точностью рисовать)
if risunok == 1
    k = 1.01:0.1:200;
    z0 = 0:0.1:400;
    s = 0:0.001:3*pi;
    sorz = 0;
    
elseif risunok == 2
    k = 0:0.001:1.5;
    z0 = 0:0.001:3;
    s = 0:0.001:5*pi;
    sorz = 0;
    
elseif risunok == 3
    k = 0:0.0001:0.25;
    z0 = 0:0.0001:0.1;
    s = 0:0.001:1*pi;
    sorz = 0;
    
elseif risunok == 4
    k = 0:0.005:1.5;
    z0 = 0:0.005:4*pi;
    s = 0:0.001:4*pi;
    sorz = 1;
    
elseif risunok == 5
    k = 0:0.001:1.5;
    z0 = 0:0.001:2;
    s = 0:0.001:1*pi;
    sorz = 1;
    
    
elseif risunok == 6
    k = 0:0.001:0.8;
    z0 = 0:0.005:2*pi;
    s = 0:0.001:1*pi;
    sorz = 1;
end

% все считается тут, вроде все нормально

[ Z0, K ] = meshgrid(z0,k);
[ res, minZ, minS, maxZ, maxS ] = result(k, z0, sorz);
[ asimpt1, asimpt2, asimpt3, asimpts, exact, sj ] = asimptotics(k, s, sorz);

res = res * 180 / pi;
minS = minS * 180 / pi;
maxS = maxS * 180 / pi;
%exact = exact * 180 / pi;

% выбор рисунка (собственно, что рисовать)

figure; 

if risunok == 1
    mesh (Z0,K,res); hold on;
    
    p = k;

    for i = 1:length(k)
       p(i) = 45; 
    end
    
    plot3(minZ,k,minS,'ro'); 
    plot3(maxZ,k,maxS,'o'); 
    temp = exact(:,3);
    plot3(temp(temp <= 400),k(temp <= 400),p(temp <= 400),'yo');
    
elseif risunok == 2

    
    [c,h]=contour (Z0,K,res,[0:5:40 41:1:44 44:0.3:46 46:1:50]); hold on;
    clabel(c,h,'manual');

    for i = 3:5
        plot(exact(:,i),k,'c','LineWidth',2); 
    end
    plot(minZ,k,'r','LineWidth',2); 
    plot(maxZ,k,'LineWidth',2); 
    [c,h]=contour (Z0,K,res,[45 45],'g','LineWidth',2);
    
elseif risunok == 3
    
    [c,h]=contour (Z0,K,res,[0:5:40 44:0.4:46]); hold on;
    clabel(c,h,'manual');
    
    [c,h]=contour (Z0,K,res,[45 45],'g','LineWidth',2);
    
    plot(exact(:,1),k,'LineWidth',2.5); hold on;
    
elseif risunok == 4
    [c,h]=contour (Z0,K,res,[0:5:40 41:1:44 44:0.3:46 46:1:50]); hold on;
    clabel(c,h,'manual');
    [c,h]=contour (Z0,K,res,[45 45],'g','LineWidth',2);

    
    for i = 1:4
    plot(exact(:,i),k,'LineWidth',2.5); hold on;
    end

    for i = 1:length(asimpts(:,1))
        plot(asimpts(i,:),k,'r','LineWidth',1.5); hold on;
    end
    for i = 1:length(sj(:,1))
        plot(sj(i,:),k,'g--','LineWidth',1.5); hold on;
    end
    
elseif risunok == 5
    plot(exact(:,1),k,'LineWidth',2.5); hold on;
    plot(asimpt1,k,'r','LineWidth',2.5); 
    plot(asimpt2,k,'g','LineWidth',2.5); 
    plot(asimpt3,k,'y','LineWidth',2.5); 

    [c,h]=contour (Z0,K,res,[0:5:40 44 46:1:50]); hold on;
    clabel(c,h,'manual');
    [c,h]=contour (Z0,K,res,[45 45],'g','LineWidth',2);

    plot(exact(:,1),k,'LineWidth',2.5); 
    plot(asimpt1,k,'r','LineWidth',2.5); 
    plot(asimpt2,k,'g','LineWidth',2.5); 
    plot(asimpt3,k,'y','LineWidth',2.5); 

    
elseif risunok == 6
    [c,h]=contour (Z0,K,res,[0:5:40 44 45.5 45.7 45:1:50]); hold on;
    clabel(c,h,'manual');
    [c,h]=contour (Z0,K,res,[45 45],'g','LineWidth',2);
    
end


% подписи

ylabel('$\kappa$','FontName','Arial Cyr','FontSize',fontsz*1.5,'Interpreter', 'latex');

if sorz == 0
    xlabel('$\eta$','FontName','Arial Cyr','FontSize',fontsz*1.5,'Interpreter', 'latex');
elseif sorz == 1
    xlabel('S','FontName','Arial Cyr','FontSize',fontsz*1.5,'Interpreter', 'latex');
end

if risunok == 1
    zlabel('Угол','FontName','Arial Cyr','FontSize',fontsz*1.5);
end

if risunok == 5
    legend({'Точное значение','1я асимпт.','2я асимпт.','3я асимпт.'},'FontSize',fontsz);
end
set (gca, 'FontSize', fontsz*1.5)

%print (gcf,'-dpng','-r150',[num2str(risunok) '.png'])
end