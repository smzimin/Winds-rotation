b = 0.01;
k = 5000;
kRight = 1;
kLeft = kRight * k;

[u,v,u1,v1, izl] = partKoefFunction(26, kRight, kLeft);

plot(u, v, 'r'); hold on;
plot([u(1,izl)-b, u(1,izl)+b], [v(1,izl)-b, v(1,izl)+b], 'r'); hold on;
plot([u(1,izl)-b, u(1,izl)+b], [v(1,izl)+b, v(1,izl)-b], 'r'); hold on;
grid;

% figure;
% 
% [u,v,u1,v1, izl] = partKoefFunction(280, kRight, kLeft);
% plot(u, v, 'b'); hold on;
% plot([u(1,izl)-b, u(1,izl)+b], [v(1,izl)-b, v(1,izl)+b], 'r'); hold on;
% plot([u(1,izl)-b, u(1,izl)+b], [v(1,izl)+b, v(1,izl)-b], 'r'); hold on;
% grid;