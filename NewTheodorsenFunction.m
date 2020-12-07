function [FF,GG] = NewTheodorsenFunction(k)
%
% Find Theodorsen's Function C(k)=F(k)+i*G(k)
%
% input: k
% output: the real and imaginary parts of C(k): FF and GG
%
% Note that because of the singularity of Bessel functions of 
% the second kind at k=0 we need to be careful when we
% calculate C(k) at k=0.
%
% Tolerance below which k is considered zero and C(k)=1.
%
epsk=0.0000000001;
%
% Value of C(k) at k=epsk
%
      J0=besselj(0,epsk);
      J1=besselj(1,epsk);
      Y0=bessely(0,epsk);
      Y1=bessely(1,epsk);
      D=(J1+Y0)^2+(Y1-J0)^2 ;                                           
      FFepsk=(J1*(J1+Y0)+Y1*(Y1-J0))/D ;
      GGepsk=-(Y1*Y0+J1*J0)/D ;
%
% Now, proceed to evaluating C(k) at the input reduced frequency k
%
if (k<=epsk)
% for very small k values (below k=epsk)
% Interpolate between the values of C(k) at k=0
% where FF0=1. and GG0=0.
% and the values at k=epsk
% (FF-FF0)/(FFepsk-FF0)=(k-0.)/(epsk-0.)
%
      FF0=1.;
      GG0=0.;
      FF=FF0+(FFepsk-FF0)*k/epsk;
      GG=GG0+(GGepsk-GG0)*k/epsk;
      %
else
      J0=besselj(0,k);
      J1=besselj(1,k);
      Y0=bessely(0,k);
      Y1=bessely(1,k);
      D=(J1+Y0)^2+(Y1-J0)^2;                                           
      FF=(J1*(J1+Y0)+Y1*(Y1-J0))/D;
      GG=-(Y1*Y0+J1*J0)/D;
end
%
% The following approximation is often used for C(k)
% It is commented out here. You may want to activate it and check
% the accuracy of the approximated C(k)
%
%      Ck=0.5+0.0075/(i*k+0.0455)+0.10055/(i*k+0.3);
%      FF=real(Ck);
%      GG=imag(Ck);
%
clear J0;
clear J1;
clear Y0;
clear Y1;
clear D;