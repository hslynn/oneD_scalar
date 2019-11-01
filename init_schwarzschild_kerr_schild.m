Globals1D;
GlobalsGR;

%metric and associated things
g00 = -(1-2.*M./x);
g01 = 2.*M./x;
g11 = 1+2.*M./x;
Pi00 = -2.*M./(x+2.*M).*sqrt(1+2.*M./x).*2.*M./power(x,2);
Pi01 = -2.*M./(x+2.*M).*sqrt(1+2.*M./x).*2.*M./power(x,2);
Pi11 = -2.*M./(x+2.*M).*sqrt(1+2.*M./x).*2.*M./power(x,2);
Phi00 = -2.*M./power(x,2);
Phi01 = -2.*M./power(x,2);
Phi11 = -2.*M./power(x,2);

%H0 = 2./power(x,2);
%H1 = 2.*(1+0)./power(x,2) + 2./x;

%deriH00 = 0.*x;
%deriH01 = 0.*x;
%deriH10 = -4./power(x,3);
%deriH11 = -2.*(0+2)./power(x,3) - 2./(x.*x);

%inverse metric
invg00 = g11./(g00.*g11-g01.*g01);
invg01 = -g01./(g00.*g11-g01.*g01);
invg11 = g00./(g00.*g11-g01.*g01);

%some auxiliary variabls
lapse = 1.0./power(-invg00, 0.5);
shift = -invg01./invg00;
normal0 = 1.0./lapse;
normal1 = -shift./lapse;
gamma11 =  1.0./g11;

%connections
gamma000 = 0.5.*shift.*Phi00 - 0.5.*lapse.*Pi00
gamma001 = 0.5.*Phi00
gamma011 = Phi01 - 0.5.*shift.*Phi11 + 0.5.*lapse.*Pi11

gamma100 = shift.*Phi01 - 0.5.*Phi00 - lapse.*Pi01
gamma101 = 0.5.*shift.*Phi11 - 0.5.*lapse.*Pi11
gamma111 = 0.5.*Phi11

gamma0 = invg00.*gamma000 + 2.*invg01.*gamma001 + invg11.*gamma011
gamma1 = invg00.*gamma100 + 2.*invg01.*gamma101 + invg11.*gamma111 - 2./x

%variables
psi =-2.*(x-x0)./power(sigma,2).*exp(-power(x-x0, 2)./power(sigma,2));
Phi_psi = (-2./power(sigma,2)+4.*power(x-x0,2)./power(sigma,4)).*exp(-power(x-x0, 2)./power(sigma,2));
Pi_psi = -normal0.*0.*x %- normal1.*Phi_psi;

