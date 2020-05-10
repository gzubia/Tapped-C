This code returns the tapped_c matching network's elements' values: 
--> Q: the required loaded Q of the resonant circuit 
--> Xp: either the inductive or capacitive reactance. They are equal at resonance 
--> Rp: the equivalent shut resistance of the inductor (Ω) 
--> C1: Tapped-C transformer's C1 capacitor (F) 
--> C2: Tapped-C transformer's C2 capacitor (F) 
--> L: Tapped-C transformer's L inductance (H) 
--> IL_dB: insertion loss of the network (dB)

For getting the above values, it requires the following inputs: 
--> fc: center frequency of the wanted resonant circuit (Hz) 
--> B: bandwidth of the wanted resonant circuit (Hz) 
--> Rs: source resistance 
--> RL: load resistance 
--> Qp: the Q of the inductor

So, the function call is,

[Q,Xp,Rp,C1,C2,L,IL_dB] = tapped_c(fc, B, Rs, Rl, Qp)

It displays the results in the command window, and it saves them a txt in the same folder it has been executed.

--------------------------------------------------------------------------------- 
Here is an example of a FM band mixer's tapped c transformer: 
---------------------------------------------------------------------------------

>> fc = 98e6; % Hz 
>> B = 20e6; % Hz 
>> Rs = 50: % Ω 
>> RL = 1.5e3; % Ω 
>> Qp = 75;

>> [Q,Xp,Rp,C1,C2,L,IL_dB] = tapped_c(fc,B,Rs,Rl,Qp);

*************************** Tapped-C resonant matching network *************************** 
| 
|   Initial data: 
| 
|       Center frequency (fc): 98.00 MHz 
|       Bandwidth (B): 20.00 MHz 
|       Serial resistance (Rs): 50.00 Ω 
|       Load resistance (Rl): 1.50 kΩ 
|       Inductor Q at 98.00 MHz (Qp): 75 
| 
| 
|   The Tapped-C transformer's characteristic equations are the following: 
| 
|       Rs' = Rs (1+C1/C2)^2 (1) 
|       CT  = C1C2/(C1+C2)   (2) 
|
|   Where 
|
|     CT is the equivalent capacitance that resonates with L| 
| 
|   First we clear C1/C2 from (1) and we will use (2) in advance to get C1 and C2 
| 
|   Let the system: 
| 
|       Qp = Rp/Xp      (3) 
|       Q  = R_total/Xp (4) 
| 
|   Resolving the proposed system, 
| 
|       Xp: 143.06 Ω 
|       Rp: 10729.59 Ω 
| 
|   Knowing that: 
| 
|       C1 = (C1/C2) * C2 (5) 
|       Ct = C1||C2 (6) 
| 
|   Then, solving the system above by replacing (5) in (6), 
| 
|         ------------------ 
|       |   L : 232.34 nH   | 
|       |   C_1: 62.18 pF   | 
|       |   C_2: 13.89 pF   | 
|        ------------------ 
| 
|   Finally, we calculate the insertion loss of the network, 
| 
|       IL = 20log10(VL_resonant/VL_noResonant) (7) 
| 
|   Solving (7)... 
| 
|       IL = 20log10((ZL*Zc1*(RL + RS))/(RL*RS*ZL + RL*RS*Zc1 + RL*RS*Zc2 + RL*ZL*Zc1 + RS*ZL*Zc1 + RS*ZL*Zc2 + RL*Zc1*Zc2 + ZL*Zc1*Zc2)) 
| 
| 
|        ------------------- 
|       |   IL(dB) : 8.84   | 
|        ------------------- 
*******************************************************************************************
