clc;clear all;close all; 
D=.1016; %with meter
epsilon=0;
A=pi/4*D^2; %with m^2
Q=2000*42/(264*24*60*60); %with m^3/sec
viscosity=0.008; %with Pa.s
density=900; %with kg/m^3
Re=density*Q*D/(A*viscosity);
f=colebrook(Re,epsilon,D)
Re=5000:188.517706056318 :100000;%Range of Reynolds Number, I determine step size to find corresponding Reynolds number in part a
epsilon=D*[0 0.002 0.004 0.006 0.008]; %Surface Roughness

for i=1:length(epsilon) %counter for different surface roughness 

    for j=1:length(Re) %to read corresponding Reynolds Number step by step

       f1(i,j)=colebrook(Re(j),epsilon(i),D); %numerical secant method
       fun=@(f) 1/sqrt(f)+2*log10(epsilon(i)/D/3.7+(2.51/(Re(j)*sqrt(f)))); %colebrook equation function 
       f2(i,j)=fzero(fun,0.1); %built function
       
    end
    
     %In order to plot
    figure(1)
    plot(Re,f1(i,:))%For x axis reynolds number and for y axis secant method function (f1) solution (i:) for every iteration values
    hold on
    
end

%to label title and legend on the given plot
figure(1)
xlabel('Reynolds Number')
ylabel('Friction Factor')
legend('for e/D=0','for e/D=0.002','for e/D=0.004','for e/D=0.006','for e/D=0.008')
title('Effects of Reynolds Number and Surface Rougness on Friction Factor')
grid

%function file for finding the friction factor
function f=colebrook(Re,e,D)
%inputs:
  %Re: Reynolds number
  %epsilon:pipes surface roughness
  %D:Pipe Diameter
%output:
  %f:Friction factor
f=[0.002 0.003];%initial values that I determined for n 
F=@(f)1/sqrt(f)+2.0*log10(e/D/3.7+2.51/(Re*sqrt(f)));
n=2;%
Tol=abs((f(n)-f(n-1))/f(n));%Tolerance value in every iterations 
while Tol>1e-6 %To reach tolerance of 10^-6 it will continue
    f(n+1)=f(n)-F(f(n))*(f(n)-f(n-1))/(F(f(n))-F(f(n-1)));%secant method equation
    n=n+1;
    Tol=abs((f(n)-f(n-1))/f(n));
    display(f(n));
end
f=f(end);
end

% Arda Pamuk, Metu Ncc Aerospace Engineering 
