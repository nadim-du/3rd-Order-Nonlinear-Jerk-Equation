%Adomian_Decomposition_Method_for_nonlinear_third_order_Jerk_Equation
close all;
clear all;
clc
nnn=27;
%Defining_JerK_Equation
syms t x(i) y  e s z C w h k yex ypp yp x2(i) zp x3(i) x4(i) r zpp D E
f(t,yp)=-yp;
g(t,y,yp,ypp)=y*yp*ypp;
%Adomian_polinomial
for j=0:nnn
y1=0;
y2=0;
y3=0;
y4=0;
for i=0:30
    y1=y1+e^i*x(i);
    y2=y2+e^i*x2(i);
    y3=y3+e^i*x3(i);
    y4=y4+e^i*x4(i);
end
f1=f(t,y1);
g1=g(t,y2,y3,y4);
a1=1/factorial(j)*diff(f1,e,j);
a2=1/factorial(j)*diff(g1,e,j);
A1=subs(a1,e,0);
A2=subs(a2,e,0);
z(1)=0.2*t;
zp(1)=0.2;
zpp(1)=0;
B1=A1;
B2=A2;
for k=0:j
    B1=subs(B1,x(k),zp(k+1));
    B2=subs(B2,x2(k),z(k+1));
    B2=subs(B2,x3(k),zp(k+1));
    B2=subs(B2,x4(k),zpp(k+1));
end
    z1(j+2)=int(int(int(B1,t),t),t);
    z2(j+2)=int(int(int(B2,t),t),t);
    z(j+2)=z1(j+2)+z2(j+2);
    zp(j+2)=diff(z(j+2),t);
    zpp(i+2)=diff(zp(j+2),t); 
end
h=0;
for i=1:nnn+2
    h=h+z(i);
    r(i)=h;
end
%Calculeted_solution_compared_with_exact_solution
aa=0;
bb=20;
hh=0.125;
n=(bb-aa)/hh;
fprintf('Time       Calculated Value    Exact Value        Absolute Error\n\n')
for i=1:n
    tt=aa+i*hh;
    tt2(i)=tt;
    yy=subs(r(nnn),t,tt);
    BB=0.2;
    ga=2*sqrt(1/(4-BB^2));
    yex=(BB/ga)*sin(ga*x)+(BB/(96*ga^3))*((-9*BB^2*ga^2-48+48*ga^2)*sin(ga*x)-BB^2*sin(3*ga*x));
    yy1=subs(yex,x,tt);
     yy2(i)=double(yy);
     yy3(i)=yy1;
    er=abs(double(yy)-double(yy1));
    fprintf('%1.4f %18.13f  %18.13f  %18.9f\n',tt,double(yy),double(yy1),double(er))
end

h1=plot(tt2, yy2, 'r', tt2, yy3,'b*')
xlabel('x value')
ylabel('y value')
leg=legend('calculated value ', 'exact value')
