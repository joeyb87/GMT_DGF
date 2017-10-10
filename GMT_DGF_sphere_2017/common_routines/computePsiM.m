%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function computes the value of the elements of Psi M using Gaussian
%  quadrature to solve the integrals. Gauss's formula for arbitrary
%  interval is taken from Abramowitz&Stegun "Handbook of mathematical
%  functions", formula 25.4.30. Abscissas and weigth factors for 32 terms
%  are taken from Table 25.4 of the same book.
%
%  INPUT:
%  myL: current value of parameter l
%  ko: mu*eps*(w_o^2) + i*mu*sigma*(w_o) from dispersion relation
%  R: radius of the spherical phantom
%
%  OUTPUT:
%  PsiEvalue: element of Psi M correspondent to l = myL
%
%  Name: computePsiM
%  Author: Riccardo Lattanzi
%  Created: 16 March 2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PsiMvalue = computePsiM(myL, ko, R)

% xi's
st_i(1)=0.0483076656;
st_i(2)=0.1444719615;
st_i(3)=0.2392873622;
st_i(4)=0.3318686022;
st_i(5)=0.4213512761;
st_i(6)=0.5068999089;
st_i(7)=0.5877157572;
st_i(8)=0.6630442669;
st_i(9)=0.7321821187;
st_i(10)=0.7944837959;
st_i(11)=0.8493676137;
st_i(12)=0.8963211557;
st_i(13)=0.9349060759;
st_i(14)=0.9647622555;
st_i(15)=0.9856115115;
st_i(16)=0.9972638618;
st_i(32:-1:17)=-st_i(1:16);
% wi's
ge_i(1)=0.096540088514727;
ge_i(2)=0.095638720079274;
ge_i(3)=0.093844399080804;
ge_i(4)=0.091173878695763;
ge_i(5)=0.087652093004403;
ge_i(6)=0.083311924226946;
ge_i(7)=0.078193895787070;
ge_i(8)=0.072345794108848;
ge_i(9)=0.065822222776361;
ge_i(10)=0.058684093478535;
ge_i(11)=0.050998059262376;
ge_i(12)=0.042835898022226;
ge_i(13)=0.034273862913021;
ge_i(14)=0.025392065309262;
ge_i(15)=0.016274394730905;
ge_i(16)=0.007018610009470;
ge_i(32:-1:17)=ge_i(1:16);

% This is based on Abramowitz&Stegun's Gauss's formula 25.4.30
x_i=(R/2)*(st_i + 1);
sb_0 = spherbessJ(myL,ko*x_i);

PsiMvalue = R/2* (sum(ge_i .* (x_i.^2 .* abs(sb_0).^2)));

% END  computePsiM.m
