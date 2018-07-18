%Initializations
Bay1 = BayConfig(1);
Bay2 = BayConfig(2);
Bay3 = BayConfig(3);
Bay4 = BayConfig(4);
Bay5 = BayConfig(5);
tube_diameters = [.01, .02, .03, .04];

D1 = tube_diameters(tube_id(1) +1); D1 = D1*2.54;
D2 = tube_diameters(tube_id(2) +1); D2 = D2*2.54;
D3 = tube_diameters(tube_id(3) +1); D3 = D3*2.54;
D4 = tube_diameters(tube_id(4) +1); D4 = D4*2.54;
D5 = tube_diameters(tube_id(5) +1); D5 = D5*2.54;

D = [D1, D2, D3, D4, D5];
length_reactor1 = 188;%cm
length_reactor2 = 99; %cm
length_bypass = 10;%cm

V = [Vr1, Vr2, Vr3, Vr4, Vr5];

%assign reactor volumes
for i =1:5
   if  new_vol_react(i) || V(i) == 0
       if new_len_react(i)
           V(i) = tube_length(i)*D(i)^2/4*3.14159*1000;
       else
           if ismember(BayConfig(i),[1,3,5])
               V(i) = length_reactor1*D(i)^2/4*3.14159*1000;
           elseif BayConfig(i) == 0
               V(i) = length_bypass*D(i)^2/4*3.14159*1000;
           elseif BayConfig(i) == 2
               V(i) = length_reactor2*D(i)^2/4*3.14159*1000;
           else
               V(i) = 1000;
           end
       end
   end
end

Vr1 = V(1);
Vr2 = V(2);
Vr3 = V(3);
Vr4 = V(4);
Vr5 = V(5);