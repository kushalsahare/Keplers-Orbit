%HW#7 problem b-Using Fourth order Runge-Kutta's method.

clc;
clear all;
t1=0;t2=pi;
N=2^15; % no of steps
h=(t2-t1)/N; %step size
t=t1:h:t2;
y=zeros(4,1,N+1);
y(:,1,1)=[0.6;0.8;-1;0];    
    
    for i=1:N
      y1=y(:,1,i);       
    y(:,1,i+1)=y(:,1,i)+h*Kepler_func(y1);
    end

plot(y(1,:),y(2,:),'--r');
hold on;
% 
% %---------------------plot ellipse------------
% 
% a=1; %horizontal radius
% b=0.8*a; %vertical radius
% x0=0.6; % x0,y0 ellipse centre coordinates
% y0=0;
% x_1=x0+a*cos(t);
% y_1=y0+b*sin(t);
% plot(x_1,y_1,'-b');
% 
% hold on;
% plot(0.6,0.8,'*r');
% title('Trajectory of comet');
% xlabel('x(\theta)');
% ylabel('y(\theta)');
% legend('Euler Method','Actual Trajectory','Initial point')



