%HW#9 problem b-Using Fourth order Runge-Kutta's Cash Karp method.

clear all;
clc;
x0=0;xf=2*pi;
TOL=5e-4;
h_max=0.25;
h_min=1e-5;
h=h_max;
x(1)=x0;
i=1;
y(:,1,1)=[0.6;0.8;1;0]; 
count=0;
    
while(x0<xf)
    
        x1=(x(i));
        y1=y(:,1,i);  
        k_1=h*Kepler_func(y1);
        %--------------------------------
        x2=(x(i)+h/5);
        y2=(y(:,1,i)+k_1/5);
        
        k_2=h*Kepler_func(y2);
        %--------------------------------
        
        x3=(x(i)+3*h/10);
        y3=(y(:,1,i)+(3*k_1/40)+(9*k_2/40));
        
        k_3=h*Kepler_func(y3);
        %-----------------------------------------------
        x4=x(i)+3*h/5;
        y4=(y(:,1,i)+ 3*k_1/10 - 9*k_2/10+ 6*k_3/5);
        
        k_4=h*Kepler_func(y4);
        
        %-----------------------------------------------
        y5=y(:,1,i)-11*k_1/5-5*k_2/2-70*k_3/27-35*k_4/27;
        x5=x(i)+h;
        k_5=h*Kepler_func(y5);
        
        %------------------------------------------------
        y6=y(:,1,i)-(1631*k_1/55296)+ (175*k_2/512) -(575*k_3/13824)+(44275*k_4/110592)-(253*k_5/4096);
        x6=x(i)+7*h/8;
        k_6=h*Kepler_func(y6);
        %--------------------------------------------------
       
        count=count+6;
        
        R = max ( abs ( -(277*k_1)/64512+(6925*k_3)/370944 -(6925*k_4)/202752 -(277*k_5)/14336+(277*k_6)/7084 ) / h );
        
        q = 0.84*( TOL / R ) ^ (1/4);
        if (R<TOL)
        w=y(:,1,i)+(2825*k_1/27648+0*k_2+ 18575*k_3/48384+13525*k_4/55296+277*k_5/14336+k_6/4); %fifth order accurate
        %w=y+(37*k_1/378+ 0*k_2+250*k_3/621+125*k_4/594+0*k_5+512*k_6/1771); % fourth order accurate
        x0=x0+h;
        i = i + 1;
        x(i) = x0;
        y(:,1,i)=w;
        end 
        h = min ( max ( q, 0.01 ), 4.0 ) * h;
	    if ( h > h_max ) 
            h = h_max;
        end;
	    if ( x0 + h > xf )
	    h = xf - x0;
	    elseif ( h < h_min )
	    disp ( 'Solution requires step size smaller than minimum' );
	    return;
	end;
        
end 

% plot(y(1,:),y(2,:),'--k');
% hold on;

%---------------------plot ellipse------------

% a=1; %horizontal radius
% b=0.8*a; %vertical radius
% x0=0.6; % x0,y0 ellipse centre coordinates
% y0=0;
% x_1=x0+a*cos(x);
% y_1=y0+b*sin(x);
% plot(x_1,y_1,'-b');
% hold on;
% plot(0.6,0.8,'*');
% xlabel('x(\theta)');
% ylabel('y(\theta)');
% title('Plot of comet''s path');
% legend('RK4 Cash Karpe Scheme','Actual Path','Start point');



