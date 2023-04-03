function[P_Kn1_real] = PIFE_rate_generator(Data)



Data.P_Kn1_real  = repmat(0,1,Data.Num_pulses);
Data.P_Kn1_real(1:30000)=repmat(0.1,1,30000);
Data.P_Kn1_real(30001:60000)=repmat(.2,1,30000);
Data.P_Kn1_real(60001:90000)=repmat(.3,1,30000);
Data.P_Kn1_real(90001:120000)=repmat(0,1,30000);
Data.P_Kn1_real(120001:150000)=repmat(1,1,30000);
Data.P_Kn1_real(150001:180000)=repmat(.6,1,30000);
Data.P_Kn1_real(180001:210000)=repmat(.7,1,30000);
Data.P_Kn1_real(210001:240000)=repmat(.8,1,30000);
Data.P_Kn1_real(240001:270000)=repmat(.9,1,30000);
Data.P_Kn1_real(270001:300000)=repmat(1,1,30000);
Data.P_Kn1_real(300001:330000)=repmat(.9,1,30000);
Data.P_Kn1_real(330001:360000)=repmat(.8,1,30000);
Data.P_Kn1_real(360001:390000)=repmat(2,1,30000);
Data.P_Kn1_real(450001:490000)=repmat(1,1,40000);
Data.P_Kn1_real(550001:560000)=repmat(2,1,10000);
Data.P_Kn1_real(650001:680000)=repmat(0,1,30000);
Data.P_Kn1_real(710001:740000)=repmat(3,1,30000);
Data.P_Kn1_real(800001:820000)=repmat(2,1,20000);
Data.P_Kn1_real(880001:910000)=repmat(2,1,30000);
Data.P_Kn1_real(910001:end)=repmat(.6,1,length(Data.P_Kn1_real(910001:end)));




P_Kn1_real  = zeros(1,Data.Num_pulses);

num00 = 0;
num03 = floor(Data.Num_pulses*0.03);
num06 = floor(Data.Num_pulses*0.06);
num09 = floor(Data.Num_pulses*0.09);
num12 = floor(Data.Num_pulses*0.12);
num15 = floor(Data.Num_pulses*0.15);
num18 = floor(Data.Num_pulses*0.18);
num21 = floor(Data.Num_pulses*0.21);
num24 = floor(Data.Num_pulses*0.24);
num27 = floor(Data.Num_pulses*0.27);
num30 = floor(Data.Num_pulses*0.30);
num33 = floor(Data.Num_pulses*0.33);
num36 = floor(Data.Num_pulses*0.36);
num39 = floor(Data.Num_pulses*0.39);
num45 = floor(Data.Num_pulses*0.45);
num49 = floor(Data.Num_pulses*0.49);
num55 = floor(Data.Num_pulses*0.55);
num56 = floor(Data.Num_pulses*0.56);
num65 = floor(Data.Num_pulses*0.65);
num68 = floor(Data.Num_pulses*0.68);
num71 = floor(Data.Num_pulses*0.71);
num74 = floor(Data.Num_pulses*0.74);
num80 = floor(Data.Num_pulses*0.80);
num82 = floor(Data.Num_pulses*0.82);
num88 = floor(Data.Num_pulses*0.88);
num91 = floor(Data.Num_pulses*0.91);
num100= floor(Data.Num_pulses*1   );

P_Kn1_real(num00+1:num03 )=repmat(0.1   ,1,num03-num00);
P_Kn1_real(num03+1:num06 )=repmat(0.2   ,1,num06-num03);
P_Kn1_real(num06+1:num09 )=repmat(0.3   ,1,num09-num06);
P_Kn1_real(num09+1:num12 )=repmat(0.0   ,1,num12-num09);
P_Kn1_real(num12+1:num15 )=repmat(1     ,1,num15-num12);
P_Kn1_real(num15+1:num18 )=repmat(0.6   ,1,num18-num15);
P_Kn1_real(num18+1:num21 )=repmat(0.7   ,1,num21-num18);
P_Kn1_real(num21+1:num24 )=repmat(0.8   ,1,num24-num21);
P_Kn1_real(num24+1:num27 )=repmat(0.9   ,1,num27-num24);
P_Kn1_real(num27+1:num30 )=repmat(1     ,1,num30-num27);
P_Kn1_real(num30+1:num33 )=repmat(0.9   ,1,num33-num30);
P_Kn1_real(num33+1:num36 )=repmat(0.8   ,1,num36-num33);
P_Kn1_real(num36+1:num39 )=repmat(2     ,1,num39-num36);
P_Kn1_real(num45+1:num49 )=repmat(1     ,1,num49-num45);
P_Kn1_real(num55+1:num56 )=repmat(2     ,1,num56-num55);
P_Kn1_real(num65+1:num68 )=repmat(0     ,1,num68-num65);
P_Kn1_real(num71+1:num74 )=repmat(3     ,1,num74-num71);
P_Kn1_real(num80+1:num82 )=repmat(2     ,1,num82-num80);
P_Kn1_real(num88+1:num91 )=repmat(2     ,1,num91-num88);
P_Kn1_real(num91+1:num100)=repmat(0.6   ,1,num100-num91);

% P_Kn1_real  = P_Kn1_real*0;

end