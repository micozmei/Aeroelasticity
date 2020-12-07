function [AA]=NewAeroMatrix(k,a,c,b,T,FF,GG,ndof)
%
% Scanlan & Rosenbaum:
% S&R pp. 199-200
%
pi=3.14159;
% the imaginary number im=i=sqrt(-1)
im=1i;
%
      Ck=FF+im*GG;
      %
      % Retrieve T functions from Vector T()
      %
      T0=T(1);
      T1=T(2);
      T2=T(3);
      T3=T(4);
      T4=T(5);
      T5=T(6);
      T6=T(7);
      T7=T(8);
      T8=T(9);
      T9=T(10);
      T10=T(11);
      T11=T(12);
      T12=T(13);
      T13=T(14);
      T14=T(15);
      T15=T(16);
      T16=T(17);
      T17=T(18);
%
% clear AA
%
AA=zeros(ndof,ndof);
%
% Scanlan & Rosenbaum, pp. 199-200
% Note The AA terms are calculated first in a form that
% makes it necessary later in function file to multiply by b and 2*pi
% Since h is distance, while alpha and beta are angles, and
% also P is force while Ma and Hb are moments, some derivations
% of the aero matrix use h/b instead of h and P*b instead of P
% to bring all entities to the same dimension
% We use such a formulation here, but correct later to return to h and P.
%
AA(1,1)=k*k-2.*im*k*Ck;
AA(1,2)=-a*k*k+(2.*(a-0.5)*Ck-1.d0)*im*k-2.*Ck;
AA(1,3)=-k*k*T1/pi+im*k*(T4/pi-T11*Ck/pi)-2.*T10/pi*Ck;
%
AA(2,1)=-k*k*a+im*k*Ck*2.*(0.5d0+a);
AA(2,2)=k*k*(a*a+0.125)+(a-0.5+2.*(0.25-a*a)*Ck)*im*k...
      + 2.*(0.5+a)*Ck;
     %
AA(2,3)=-T15/pi-im*k*T16/pi+2.*k*k*T13/pi+ ...
     2.*(a+0.5)*Ck*T10/pi+2.*(a+0.5)*Ck*im*k*T11/2./pi;
%
AA(3,1)=-k*k*(T1/pi)-im*k*(T12/pi)*Ck;
%
AA(3,3)=-k*k*T3/pi/pi+im*k*...
   (T4*T11/2./pi/pi -T11*T12/2./pi/pi*Ck )...
    -(T5-T4*T10)/pi/pi -T10*T12/pi/pi*Ck;
%
AA(3,2)=-im*k*T17/1./pi+2.*T13*k*k/pi-2.*Ck*T12/2./pi...
         -Ck*T12/pi*im*k*(0.5-a);
%
% Now multiply the whole AA matrix above by 2pi
%
AA=2.*pi*b^2*AA;
%
% and then:
%
% Modify to account for P*b --> P and h/b --> h 
%
for ii=1:ndof
    AA(1,ii)=AA(1,ii)/b;
    AA(ii,1)=AA(ii,1)/b;
end   