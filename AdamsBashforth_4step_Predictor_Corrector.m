%HW#9 problem a-Using Adams-bashforth scheme

clear all;
clc;
t1=0;t2=2*pi;
N=2^10; % no of steps
h=(t2-t1)/N; %step size
t=t1:h:t2;
y=zeros(4,1,N+1);
y(:,1,1)=[0.6;0.8;1;0];    
    
    for i=1:3
        %----------------------------------
        x1=t(i);
        y1=y(:,1,i);
        k_1=h*Kepler_func(y1);
        %----------------------------------
        x2=t(i)+h/2;
        y2=y(:,1,i)+k_1/2;
        k_2=h*Kepler_func(y2);
        %----------------------------------
        x3=t(i)+h/2;
        y3=y(:,1,i)+k_2/2;
        k_3=h*Kepler_func(y3);
        %---------------------------------
        x4=t(i)+h;
        y4=y(:,1,i)+k_3;
        k_4=h*Kepler_func(y4);
       
        y(:,1,i+1)=y(:,1,i)+(k_1+2*k_2+2*k_3+k_4)/6;
    end
    
    for i=4:N
        w=y(:,1,i)+h*(55*Kepler_func(y(:,1,i))/24-59*Kepler_func(y(:,1,i))/24+37*Kepler_func(y(:,1,i))/24-9*Kepler_func(y(:,1,i-3))/24);
        y(:,1,i+1)=y(:,1,i)+h*((9*Kepler_func(w)+19*Kepler_func(y(:,1,i))-5*Kepler_func(y(:,1,i))+Kepler_func(y(:,1,i-3)))/24);
    end
plot(y(1,:),y(2,:),'--k');
hold on;

%---------------------plot ellipse------------

a=1; %horizontal radius
b=0.8*a; %vertical radius
x0=0.6; % x0,y0 ellipse centre coordinates
y0=0;
x_1=x0+a*cos(t);
y_1=y0+b*sin(t);
plot(x_1,y_1,'-b');

hold on;
plot(0.6,0.8,'*');
title('Plot of comet''s path');
xlabel('x(\theta)');
ylabel('y(\theta)');
legend('4 step Predictor-Corrector Schemes','Actual path','Start point');



