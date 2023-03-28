function[F_Kn1_real] = FRET_rate_generator(Data)


F_Kn1_real  = zeros(1,Data.Num_pulses);

num00 = 0;
num03 = floor(Data.Num_pulses*0.03);
num06 = floor(Data.Num_pulses*0.06);
num09 = floor(Data.Num_pulses*0.09);
num15 = floor(Data.Num_pulses*0.15);
num16 = floor(Data.Num_pulses*0.16);
num20 = floor(Data.Num_pulses*0.20);
num23 = floor(Data.Num_pulses*0.23);
num28 = floor(Data.Num_pulses*0.28);
num35 = floor(Data.Num_pulses*0.35);
num40 = floor(Data.Num_pulses*0.40);
num50 = floor(Data.Num_pulses*0.50);
num60 = floor(Data.Num_pulses*0.60);
num65 = floor(Data.Num_pulses*0.65);
num80 = floor(Data.Num_pulses*0.80);
num82 = floor(Data.Num_pulses*0.82);
num90 = floor(Data.Num_pulses*0.90);
num93 = floor(Data.Num_pulses*0.93);

F_Kn1_real(num00+1:num03)=repmat(1.01  ,1,num03-num00);
F_Kn1_real(num06+1:num09)=repmat(2     ,1,num09-num06);
F_Kn1_real(num15+1:num16)=repmat(4     ,1,num16-num15);
F_Kn1_real(num20+1:num23)=repmat(3     ,1,num23-num20);
F_Kn1_real(num28+1:num35)=repmat(1.01  ,1,num35-num28);
F_Kn1_real(num40+1:num50)=repmat(2     ,1,num50-num40);
F_Kn1_real(num60+1:num65)=repmat(6     ,1,num65-num60);
F_Kn1_real(num80+1:num82)=repmat(3     ,1,num82-num80);
F_Kn1_real(num90+1:num93)=repmat(4     ,1,num93-num90);


end