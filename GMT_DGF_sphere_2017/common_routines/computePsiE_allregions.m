%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function computes the value of the elements of Psi E using Gaussian
%  quadrature to solve the integrals. Gauss's formula for arbitrary
%  interval is taken from Abramowitz&Stegun "Handbook of mathematical
%  functions", formula 25.4.30. Abscissas and weigth factors for 32 terms
%  are taken from Table 25.4 of the same book.
%
%  INPUT:
%  myL: current value of parameter l
%  ko: mu*eps*(w_o^2) + i*mu*sigma*(w_o) from dispersion relation
%  r1,r2: extremes of the interval of integration
%  type: 1 --> Bessel function, 2 and > 2 --> Hankel function of the first kind
%
%  OUTPUT:
%  PsiEvalue: element of Psi E correspondent to l = myL
%
%  Name: computePsiE
%  Author: Riccardo Lattanzi
%  Created: 16 March 2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PsiEvalue = computePsiE_allregions(myL, ko, r1 ,r2, type)

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

if ~( (type ==  1) || (type ==  2) || (type ==  3)|| (type ==  4) )
    type = 1;
    disp('-------------------------------------------------------------------------------------');
    disp('**WARNING** type is not 1,2,3 or 4, so the PsiE integral is calculated using type = 1');
    disp('-------------------------------------------------------------------------------------');
end

% This is based on Abramowitz&Stegun's Gauss's formula 25.4.30
x_i=((r2 - r1)/2)*st_i + (r2 + r1)/2;

switch type
    case 1
        sb_0 = spherbessJ(myL,ko*x_i);
        sb_d = derRspherbessJ(myL,ko*x_i);
        PsiEvalue = ((r2 - r1)/2)*(1/(abs(ko)^2))*( sum(ge_i.*abs(sb_d).^2) + myL*(myL+1)*sum(ge_i.*abs(sb_0).^2));
    case 2
        sb_0 = spherbessH(myL,ko*x_i);
        sb_d = derRspherbessH(myL,ko*x_i);
        PsiEvalue = ((r2 - r1)/2)*(1/(abs(ko)^2))*( sum(ge_i.*abs(sb_d).^2) + myL*(myL+1)*sum(ge_i.*abs(sb_0).^2));
    case 3
        sb_0_1 = spherbessJ(myL,ko*x_i);
        sb_0_2 = conj(spherbessH(myL,ko*x_i));
        sb_d_1 = derRspherbessJ(myL,ko*x_i);
        sb_d_2 = conj(derRspherbessH(myL,ko*x_i));
        PsiEvalue = ((r2 - r1)/2)*(1/(abs(ko)^2))*( sum(ge_i.*(sb_d_1.*sb_d_2)) + myL*(myL+1)*sum(ge_i.*(sb_0_1.*sb_0_2)));
    case 4
        sb_0_1 = conj(spherbessJ(myL,ko*x_i));
        sb_0_2 = spherbessH(myL,ko*x_i);
        sb_d_1 = conj(derRspherbessJ(myL,ko*x_i));
        sb_d_2 = derRspherbessH(myL,ko*x_i);
        PsiEvalue = ((r2 - r1)/2)*(1/(abs(ko)^2))*( sum(ge_i.*(sb_d_1.*sb_d_2)) + myL*(myL+1)*sum(ge_i.*(sb_0_1.*sb_0_2)));
end
        
% ascissa = [0.0001:0.0005:0.15];
% figure; plot(ascissa, abs(ko.*derspherbess(0,ko*ascissa)).^2);

% PsiEvalue = ((r2 - r1)/2)*(1/(abs(ko)^2))*( sum(ge_i.*abs(sb_d).^2) + myL*(myL+1)*sum(ge_i.*abs(sb_0).^2));
% PsiEvalue = (R/2)*( sum(ge_i.*abs(sb_d).^2) + (2/abs(ko))*sum(ge_i.*abs(sb_d).*abs(sb_0)) + ...
%     (1/(abs(ko)^2))*(myL*(myL+1) + 1)*sum(ge_i.*abs(sb_0).^2));

% END  computePsiE.m
