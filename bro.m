%Givna värden 
Lh=3;   Lb=4;
E=210e9;
I=24.92e-6;  %Balk
A1=5.425e-3; %Balk
A2=1.2e-3;   %Stång
q=30e3;
h=0.160;

%Matris
K=zeros(23);
f=zeros(23,1);

%koordinater
Coord=[0 0; Lb 0; 2*Lb 0; 3*Lb 0; 3.5*Lb 0; 
    0.5*Lb Lh; 1.5*Lb Lh; 2.5*Lb Lh; 3.5*Lb Lh];
Dof1=[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15];
Dof2=[1 2; 4 5; 7 8; 10 11; 13 14; 16 17; 18 19; 20 21; 22 23];

%Element
nel=2;
ep1=[E A1 I];
ep2=[E A2];

Edof1=[1 1 2 3 4 5 6;
    2 4 5 6 7 8 9;
    3 7 8 9 10 11 12;
    4 10 11 12 13 14 15];

Edof2=[5 16 17 18 19;
    6 18 19 20 21;
    7 20 21 22 23;
    8 1 2 16 17;
    9 16 17 4 5;
    10 4 5 18 19;
    11 18 19 7 8;
    12 7 8 20 21;
    13 20 21 10 11;
    14 10 11 22 23];

[Ex1,Ey1] = coordxtr(Edof1,Coord,Dof1,nel);    
[Ex2,Ey2] = coordxtr(Edof2,Coord,Dof2,nel);  
    
eldraw2(Ex1,Ey1,[1 3 1]);
eldraw2(Ex2,Ey2,[1 2 1]);

%assemblera elementmatrisen
eq=[0 -q];

for i=1:4
    [Ke,fe]=beam2e(Ex1(i,:),Ey1(i,:),ep1,eq);
    [K,f]=assem(Edof1(i,:),K,Ke,f,fe);
end

for i=1:10
    Ke=bar2e(Ex2(i,:),Ey2(i,:),ep2);
    K=assem(Edof2(i,:),K,Ke);
end

%3a
%lösa ekvationssystemet
bc=[2 0; 13 0; 15 0; 22 0];

[a,r]=solveq(K,f,bc)


Ed1=extract(Edof1,a); 
Ed2=extract(Edof2,a);   

[sfac]=scalfact2(Ex1,Ey1,Ed1);
eldisp2(Ex1,Ey1,Ed1,[2 1 1],sfac);
eldisp2(Ex2,Ey2,Ed2,[2 1 1],sfac);

%3b
%Balk
es1 = beam2s(Ex1(1,:),Ey1(1,:),ep1,Ed1(1,:),eq,39);
es2 = beam2s(Ex1(2,:),Ey1(2,:),ep1,Ed1(2,:),eq,39);
es3 = beam2s(Ex1(3,:),Ey1(3,:),ep1,Ed1(3,:),eq,39);
es4 = beam2s(Ex1(4,:),Ey1(4,:),ep1,Ed1(4,:),eq,39);

%Stång
es5 = bar2s(Ex2(1,:),Ey2(1,:),ep2,Ed2(1,:));
es6 = bar2s(Ex2(2,:),Ey2(2,:),ep2,Ed2(2,:));
es7 = bar2s(Ex2(3,:),Ey2(3,:),ep2,Ed2(3,:));
es8 = bar2s(Ex2(4,:),Ey2(4,:),ep2,Ed2(4,:));
es9 = bar2s(Ex2(5,:),Ey2(5,:),ep2,Ed2(5,:));
es10 = bar2s(Ex2(6,:),Ey2(6,:),ep2,Ed2(6,:));
es11 = bar2s(Ex2(7,:),Ey2(7,:),ep2,Ed2(7,:));
es12 = bar2s(Ex2(8,:),Ey2(8,:),ep2,Ed2(8,:));
es13 = bar2s(Ex2(9,:),Ey2(9,:),ep2,Ed2(9,:));
es14 = bar2s(Ex2(10,:),Ey2(10,:),ep2,Ed2(10,:));

%Rita moment
sfac=1e-4;
figure(2)
plotpar=[4 1];
eldia2(Ex1(1,:),Ey1(1,:),es1(:,3),plotpar,sfac);
eldia2(Ex1(2,:),Ey1(2,:),es2(:,3),plotpar,sfac);
eldia2(Ex1(3,:),Ey1(3,:),es3(:,3),plotpar,sfac);
eldia2(Ex1(4,:),Ey1(4,:),es4(:,3),plotpar,sfac);
xlabel('längd (m)')
ylabel('Moment (Nm 1e-4)')
title('Momentdiagram')

%Tvärkraft
sfac=5e-5;
figure(3)
plotpar=[3 1];
eldia2(Ex1(1,:),Ey1(1,:),es1(:,2),plotpar,sfac);
eldia2(Ex1(2,:),Ey1(2,:),es2(:,2),plotpar,sfac);
eldia2(Ex1(3,:),Ey1(3,:),es3(:,2),plotpar,sfac);
eldia2(Ex1(4,:),Ey1(4,:),es4(:,2),plotpar,sfac);
xlabel('längd (m)')
ylabel('Tvärkraft (N 5e-5)')
title('Tvärkraftsdiagram')

%Normalkraftsfigur ritas
sfac=5e-6;
figure(4)
plotpar=[2 1];
eldia2(Ex1(1,:),Ey1(1,:),es1(:,1),plotpar,sfac);
eldia2(Ex1(2,:),Ey1(2,:),es2(:,1),plotpar,sfac);
eldia2(Ex1(3,:),Ey1(3,:),es3(:,1),plotpar,sfac);
eldia2(Ex1(4,:),Ey1(4,:),es4(:,1),plotpar,sfac);
% eldia2(Ex2(1,:),Ey2(1,:),es5,plotpar,sfac);
% eldia2(Ex2(2,:),Ey2(2,:),es6,plotpar,sfac);
% eldia2(Ex2(3,:),Ey2(3,:),es7,plotpar,sfac);
% eldia2(Ex2(4,:),Ey2(4,:),es8,plotpar,sfac);
% eldia2(Ex2(5,:),Ey2(5,:),es9,plotpar,sfac);
% eldia2(Ex2(6,:),Ey2(6,:),es10,plotpar,sfac);
% eldia2(Ex2(7,:),Ey2(7,:),es11,plotpar,sfac);
% eldia2(Ex2(8,:),Ey2(8,:),es12,plotpar,sfac);
% eldia2(Ex2(9,:),Ey2(9,:),es13,plotpar,sfac);
% eldia2(Ex2(10,:),Ey2(10,:),es14,plotpar,sfac);
xlabel('längd (m)')
ylabel('Normalkraft (N 5e-6)')
title('Normalkraftsdiagram')


%4
%Givna värden 
A4=1e-20;   %uppg4 försumbararea

%Matris
K=zeros(23);
f=zeros(23,1);

%koordinater
Coord=[0 0; Lb 0; 2*Lb 0; 3*Lb 0; 3.5*Lb 0; 
    0.5*Lb Lh; 1.5*Lb Lh; 2.5*Lb Lh; 3.5*Lb Lh];
Dof1=[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15];
Dof2=[1 2; 4 5; 7 8; 10 11; 13 14; 16 17; 18 19; 20 21; 22 23];

%Element
nel=2;
ep1=[E A1 I];
ep4=[E A4];

Edof1=[1 1 2 3 4 5 6;
    2 4 5 6 7 8 9;
    3 7 8 9 10 11 12;
    4 10 11 12 13 14 15];

Edof2=[5 16 17 18 19;
    6 18 19 20 21;
    7 20 21 22 23;
    8 1 2 16 17;
    9 16 17 4 5;
    10 4 5 18 19;
    11 18 19 7 8;
    12 7 8 20 21;
    13 20 21 10 11;
    14 10 11 22 23];

[Ex1,Ey1] = coordxtr(Edof1,Coord,Dof1,nel);    
[Ex2,Ey2] = coordxtr(Edof2,Coord,Dof2,nel);  

figure(5)
eldraw2(Ex1,Ey1,[1 3 1]);
eldraw2(Ex2,Ey2,[1 2 1]);

%assemblera elementmatrisen
eq=[0 -q];

for i=1:4
    [Ke,fe]=beam2e(Ex1(i,:),Ey1(i,:),ep1,eq);
    [K,f]=assem(Edof1(i,:),K,Ke,f,fe);
end

for i=1:10
    Ke=bar2e(Ex2(i,:),Ey2(i,:),ep4);
    K=assem(Edof2(i,:),K,Ke);
end

bc=[2 0; 13 0; 15 0; 22 0];

[a4,r4]=solveq(K,f,bc);

Ed41=extract(Edof1,a4); 
Ed42=extract(Edof2,a4);   

[sfac]=scalfact2(Ex1,Ey1,Ed41);
eldisp2(Ex1,Ey1,Ed41,[2 1 1],sfac);
eldisp2(Ex2,Ey2,Ed42,[2 1 1],sfac);

%största utböjningen
vmax=max(abs(Ed41))


%5
%Max krafter och jämnvikt
esbar=[es5; es6; es7; es8 ; es9; es10; es11; es12; es13; es14];
esbeam=[es1; es2; es3; es4];

%5a
A2=1.2e-3;
[maxNormalkraft_5a, Stang_5a] = max(abs(esbar(:))) 
%max normalkraft i stängerna
%stang är vilken stång (där element motsvaras av stång+4)
maxNormalspanning_5a=maxNormalkraft_5a/A2

%5b
[maxM, Mx] = max(abs(esbeam(:,3))); %maxM=största moment, Mx=plats maxM inträffar
maxMoment_5b=esbeam(Mx,3)
maxMx_5b=Mx*(Lb/Mx) %x-koordinat för Max moment

%5c
maxN = esbeam(Mx,1);
y=linspace(-h/2,h/2,10);
Spanningsfordelning=(maxN/A1)+((maxM/I)*y)
figure(6)
plot(Spanningsfordelning,y)
xline(0,'--r');
title('spänningsfördelning 5c')


%6, Deformation av tredje balkelement
x=linspace(0,Lb,39);
vp=(-q/(E*I)*(((x.^4)/24)-((Lb*x.^3)/12)+(((Lb.^2)*x.^2)/24)));
N=[1-(3*((x.^2)/(Lb.^2)))+(2*((x.^3)/(Lb.^3))) 
    x-(2*((x.^2)/Lb))+((x.^3)/(Lb.^2))
    3*((x.^2)/(Lb.^2))-2*((x.^3)/(Lb.^3))
    -((x.^2)/Lb)+((x.^3)/(Lb.^2))];
a3=[a(8);a(9);a(11);a(12)];
v=N'*a3+vp';
figure(7)
plot(x,v)
grid on
title('utböjning för balkelement 3 [m]')


%8, Icke-linjäritet
L=Lb;
%Matris
K=zeros(15);
f=zeros(15,1);

%koordinater
Coord=[0 0; L 0; 2*L 0; 3*L 0; 3.5*L 0; 
    0.5*L Lh; 1.5*L Lh; 2.5*L Lh; 3.5*L Lh];
Dof1=[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15];
Dof2=[1 2; 4 5; 7 8; 10 11; 13 14; 16 17; 18 19; 20 21; 22 23];

%Element
nel=2;
ep1=[E A1 I];
ep2=[E A2];

Edof1=[1 1 2 3 4 5 6;
    2 4 5 6 7 8 9;
    3 7 8 9 10 11 12;
    4 10 11 12 13 14 15];

Edof2=[5 16 17 18 19;
    6 18 19 20 21;
    7 20 21 22 23;
    8 1 2 16 17;
    9 16 17 4 5;
    10 4 5 18 19;
    11 18 19 7 8;
    12 7 8 20 21;
    13 20 21 10 11;
    14 10 11 22 23];

[Ex1,Ey1] = coordxtr(Edof1,Coord,Dof1,nel);    
[Ex2,Ey2] = coordxtr(Edof2,Coord,Dof2,nel);  

figure(8)
eldraw2(Ex1,Ey1,[1 3 1]);
eldraw2(Ex2,Ey2,[1 2 1]);

%assemblera elementmatrisen
eq=[0 -q];

eps = 0.0001;
N = [0.01 0 0 0];
N0 = [ 1 1 1 1];
n=0;

while (abs((N(1)-N0(1))/N0(1)) > eps)
    n = n+1;

K = zeros(23);
f1 = zeros(23,1);

for i = 1:4
[Ke1, fe1] = beam2g(Ex1(i,:),Ey1(i,:),ep1,N(i),-q);
[K,f1] = assem(Edof1(i,:),K,Ke1,f1,fe1);
end

for i = 1:10
    Ke2= bar2e(Ex2(i,:),Ey2(i,:),ep2);
    K = assem(Edof2(i,:),K,Ke2);
end

% Randvillkor:

bc = [1 (0.075/2);
    2 0;
    13 0;
    15 0;
    22 0];

% Lös ekvationen:

[a,r] = solveq(K,f1,bc);
Ed1 = extract(Edof1,a);
Ed2 = extract(Edof2,a);

Es1 = beam2gs(Ex1(1,:),Ey1(1,:),ep1,Ed1(1,:),N(1),-q);
Es2 = beam2gs(Ex1(2,:),Ey1(2,:),ep1,Ed1(2,:),N(2),-q);
Es3 = beam2gs(Ex1(3,:),Ey1(3,:),ep1,Ed1(3,:),N(3),-q);
Es4 = beam2gs(Ex1(4,:),Ey1(4,:),ep1,Ed1(4,:),N(4),-q);

N0 = N;
N = [Es1(1,1) Es2(1,1) Es3(1,1) Es4(1,1)];
if (n>20)
        disp('the solution doesn`t converge')
        return
end

disp("Iterations number")
disp(n)
end

[sfac] = scalfact2(Ex1,Ey1,Ed1, 0.1);
eldisp2(Ex1,Ey1,Ed1,[3 1 1],sfac);
eldisp2(Ex2,Ey2,Ed2,[3 1 1],sfac)

esbeamI=[Es1; Es2; Es3; Es4];
[maxM, Mx] = max(abs(esbeamI(:,3))); %maxM=största moment, Mx=plats maxM inträffar
maxMoment_8=maxM


%9, Fjädrande upplag [Eftergivlig]
%Givna värden 
Lh=3;   Lb=4;
E=210e9;
I=24.92e-6;  %Balk
A1=5.425e-3; %Balk
A2=1.2e-3;   %Stång
q=30e3;
%Givna värden för fundament
kx=10e6;    %axiell fundamentstyvhet
ky=10e6;    %transversel fundamentstyvhet
A3=1;       %fundamentets area
hf=1;       %fundamentets höjd
If=A3/12;   %Tröghetsmoment för fundament

%Matris
K=zeros(23);
f=zeros(23,1);

%koordinater
Coord=[0 0; Lb 0; 2*Lb 0; 3*Lb 0; 3.5*Lb 0; 
    0.5*Lb Lh; 1.5*Lb Lh; 2.5*Lb Lh; 3.5*Lb Lh];
Dof1=[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15];
Dof2=[1 2; 4 5; 7 8; 10 11; 13 14; 16 17; 18 19; 20 21; 22 23];

%Element
nel=2;
ep1=[E A1 I];
ep2=[E A2];

Edof1=[1 1 2 3 4 5 6;
    2 4 5 6 7 8 9;
    3 7 8 9 10 11 12;
    4 10 11 12 13 14 15];

Edof2=[5 16 17 18 19;
    6 18 19 20 21;
    7 20 21 22 23;
    8 1 2 16 17;
    9 16 17 4 5;
    10 4 5 18 19;
    11 18 19 7 8;
    12 7 8 20 21;
    13 20 21 10 11;
    14 10 11 22 23];

Edof3=[15 1 2 3];

[Ex1,Ey1] = coordxtr(Edof1,Coord,Dof1,nel);    
[Ex2,Ey2] = coordxtr(Edof2,Coord,Dof2,nel);  

figure(10)
eldraw2(Ex1,Ey1,[1 3 1]);
eldraw2(Ex2,Ey2,[1 2 1]);
axis([0 Lb -1 Lh])

%assemblera elementmatrisen
eq=[0 -q];

for i=1:4
    [Ke,fe]=beam2e(Ex1(i,:),Ey1(i,:),ep1,eq);
    [K,f]=assem(Edof1(i,:),K,Ke,f,fe);
end

%assemblera in fundamentet
Ke3=[kx*A3 0 kx*A3*hf;
    0 ky*A3 0;
    kx*A3*hf 0 ky*If+kx*A3*h.^2];
K=assem(Edof3(1,:),K,Ke3);

for i=1:10
    Ke=bar2e(Ex2(i,:),Ey2(i,:),ep2);
    K=assem(Edof2(i,:),K,Ke);
end

%lösa ekvationssystemet
bc=[13 0; 15 0; 22 0];

[a,r]=solveq(K,f,bc);

Ed1=extract(Edof1,a); 
Ed2=extract(Edof2,a);   

[sfac]=scalfact2(Ex1-1,Ey1,Ed1);
eldisp2(Ex1,Ey1,Ed1,[2 4 1],sfac);
eldisp2(Ex2,Ey2,Ed2,[2 4 1],sfac);

%9, Fjädrande upplag [Oeftergivlig]
%Matris
K=zeros(23);
f=zeros(23,1);

%koordinater
Coord=[0 0; Lb 0; 2*Lb 0; 3*Lb 0; 3.5*Lb 0; 
    0.5*Lb Lh; 1.5*Lb Lh; 2.5*Lb Lh; 3.5*Lb Lh];
Dof1=[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15];
Dof2=[1 2; 4 5; 7 8; 10 11; 13 14; 16 17; 18 19; 20 21; 22 23];

%Element
nel=2;
ep1=[E A1 I];
ep2=[E A2];

Edof1=[1 1 2 3 4 5 6;
    2 4 5 6 7 8 9;
    3 7 8 9 10 11 12;
    4 10 11 12 13 14 15];

Edof2=[5 16 17 18 19;
    6 18 19 20 21;
    7 20 21 22 23;
    8 1 2 16 17;
    9 16 17 4 5;
    10 4 5 18 19;
    11 18 19 7 8;
    12 7 8 20 21;
    13 20 21 10 11;
    14 10 11 22 23];

[Ex1,Ey1] = coordxtr(Edof1,Coord,Dof1,nel);    
[Ex2,Ey2] = coordxtr(Edof2,Coord,Dof2,nel);  

%assemblera elementmatrisen
eq=[0 -q];

for i=1:4
    [Ke,fe]=beam2e(Ex1(i,:),Ey1(i,:),ep1,eq);
    [K,f]=assem(Edof1(i,:),K,Ke,f,fe);
end

for i=1:10
    Ke=bar2e(Ex2(i,:),Ey2(i,:),ep2);
    K=assem(Edof2(i,:),K,Ke);
end

%lösa ekvationssystemet
bc=[1 0; 2 0; 3 0; 13 0; 15 0; 22 0];

[a,r]=solveq(K,f,bc);

Ed3=extract(Edof1,a); 
Ed4=extract(Edof2,a);   

[sfac]=scalfact2(Ex1,Ey1,Ed1);
eldisp2(Ex1,Ey1,Ed3,[2 1 1],sfac);
eldisp2(Ex2,Ey2,Ed4,[2 1 1],sfac);

