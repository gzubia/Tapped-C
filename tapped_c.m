function [Q,Xp,Rp,C1,C2,L,IL_dB] = tapped_c(fc,B,Rs,Rl,Qp)
% This code returns the tapped_c matching network's elements' values:
                                    % Q:     the required loaded Q of the resonant circuit
                                    % Xp:    either the inductive or capacitive reactance. They are equal at resonance (Ω)
                                    % Rp:    the equivalent shut resistance of the inductor (Ω)
                                    % C1:    Tapped-C transformer's C1 capacitor (F)
                                    % C2:    Tapped-C transformer's C2 capacitor (F)
                                    % L:     Tapped-C transformer's L inductance (H)
                                    % IL_dB: the network's insertion loss (dB)
                                    
%%% Code created by Gorka Zubia Garea the 10-06-2020
name = sprintf('Tapped-C_fc%4.2f_B%4.2f_Rs%4.2f_Rl%4.2f_Qp%d.txt',fc,B,Rs,Rl,Qp);
fid=fopen(name,'w'); %Create txt

%% Ploteatu saraearen hasierako datuak
fprintf(fid,'\n*************************** Tapped-C resonant matching network ***************************\n');
fprintf('\n*************************** Tapped-C resonant matching network ***************************\n');
[fc_p,fc_u]  = units(fc); % Hz
[B_p,  B_u]  = units(B);  % Hz
[Rs_p, Rs_u] = units(Rs); % Ω
[Rl_p, Rl_u] = units(Rl); % Ω
fprintf('|\n|\tInitial data:\n|\n|\t\tCenter frequency (fc): %12.2f %cHz\n|\t\tBandwidth (B): %20.2f %cHz\n|\t\tSerial resistance (Rs): %11.2f %cΩ\n|\t\tLoad   resistance (Rl): %11.2f %cΩ\n|\t\tInductor Q at %4.2f %cHz (Qp): %d\n|',fc_p,fc_u,B_p,B_u,Rs_p,Rs_u,Rl_p,Rl_u,fc_p,fc_u,Qp);
fprintf(fid,'|\n|\tInitial data:\n|\n|\t\tCenter frequency (fc): %12.2f %cHz\n|\t\tBandwidth (B): %20.2f %cHz\n|\t\tSerial resistance (Rs): %11.2f %cΩ\n|\t\tLoad   resistance (Rl): %11.2f %cΩ\n|\t\tInductor Q at %4.2f %cHz (Qp): %d\n|',fc_p,fc_u,B_p,B_u,Rs_p,Rs_u,Rl_p,Rl_u,fc_p,fc_u,Qp);
%% Calculate the baseline
%%% We know that,
%%%                 Rs' = Rs (1+C1/C2)^2 (1)
%%%
%%%                 CT  = C1C2/(C1+C2)   (2)
fprintf("\n|\n|\tThe Tapped-C transformer's characteristic equations are the following:\n|\n|");
fprintf(fid,"\n|\n|\tThe Tapped-C transformer's characteristic equations are the following:\n|\n|");

fprintf("\t\t\tRs' = Rs (1+C1/C2)^2 (1)\n|\t\t\tCT  = C1C2/(C1+C2)   (2)\n|");
fprintf(fid,"\t\t\tRs' = Rs (1+C1/C2)^2 (1)\n|\t\t\tCT  = C1C2/(C1+C2)   (2)\n|");

fprintf("\tWhere\n|\t\tCT is the equivalent capacitance that resonates with L");
fprintf(fid,"\tWhere\n|\t\tCT is the equivalent capacitance that resonates with L");

%%% So...
Rs_prima = Rl;                    % Ohm
c1_c2    = sqrt(Rs_prima/Rs) - 1; % c1/c2
w        = 2*pi*fc;               % radHz
fprintf("|\n|\n|\tFirst we clear C1/C2 from (1) and we will use (2) in advance to get C1 and C2\n");
fprintf(fid,"|\n|\n|\tFirst we clear C1/C2 from (1) and we will use (2) in advance to get C1 and C2\n");

%% Calculate the Tapped-C resonant matching network
%%%%%%%%%%%%%%%%%%%%%%
%      Loaded Q      %
%%%%%%%%%%%%%%%%%%%%%%
Q = fc/B;
%%%%%%%%%%%%%%%%%%%%%%
%      Xp and Rp     %
%%%%%%%%%%%%%%%%%%%%%%
%%% Get Xp eta Rp values
%%%% Rtotal = Rsprima||Rl||Rp
syms Rp
R_total = 1/((1/Rs_prima)+(1/Rp)+(1/Rl));
%%%% Qp = Rp/Qp
syms Xp
eqn_1 = Qp == Rp/Xp;
%%%% Q  = Rtotal/Xp
eqn_2 = Q == R_total/Xp;
%%%% Solve the system
S   = solve([eqn_1, eqn_2],[Rp Xp]);
Xp  = S.Xp;
Rp  = S.Rp;
%%%% Display the Xp eta Rp values
fprintf('|\n|\tLet the system:\n|\n|\t\t\tQp = Rp/Xp\t (3)\n|\t\t\tQ  = R_total/Xp\t (4)\n');
fprintf(fid,'|\n|\tLet the system:\n|\n|\t\t\tQp = Rp/Xp\t (3)\n|\t\t\tQ  = R_total/Xp\t (4)\n');

fprintf('|\n|\tResolving the proposed system,\n|\n|\t\t\tXp:%9.2f Ω\n|\t\t\tRp:%9.2f Ω\n|\n|',Xp,Rp);
fprintf(fid,'|\n|\tResolving the proposed system,\n|\n|\t\t\tXp:%9.2f Ω\n|\t\t\tRp:%9.2f Ω\n|\n|',Xp,Rp);

%%%%%%%%%%%%%%%%%%%%%%
%          L         %
%%%%%%%%%%%%%%%%%%%%%%
%%% Get L
%%%% Knowing that Xp = wL = 2πfcL
L            = Xp/w;     % H. Return L in the original data format
[L_adap,L_u] = units(L); % Adapt the returned L to show in convinient unit measure (i.e nH, mH, ...)

%%%%%%%%%%%%%%%%%%%%%%
%      C1 and C2     %
%%%%%%%%%%%%%%%%%%%%%%
%%% Get C1 and C2
%%%% Assuming that, Xp = 1/(wCt), where Ct = C1||C2
C_t = 1/(w*Xp);                         % (F). Get total C
syms C_1 C_2
eqn_C1 = C_1 == c1_c2 * C_2;            % C1 is C1/C2 * C2
eqn_C2 = C_t == 1/sum((1/C_1)+(1/C_2)); % Total C is the series equivalent capacitance of C1 and C2
%%%%% Solving the proposed equation
S_c = solve([eqn_C1,eqn_C2],[C_1,C_2]);
C1 = S_c.C_1;                           % (F). C1 value in original units
C2 = S_c.C_2;                           % (F). C2 value in original units
[C_1,C_1_u] = units(S_c.C_1);           % C1 value in convinient units (i.e. pF, nF, ...)
[C_2,C_2_u] = units(S_c.C_2);           % C2 value in convinient units (i.e. pF, nF, ...)
fprintf('\tKnowing that:\n|\n|\t\t\tC1 = (C1/C2) * C2   (5)\n|\t\t\tCt = C1||C2\t    (6)\n');
fprintf(fid,'\tKnowing that:\n|\n|\t\t\tC1 = (C1/C2) * C2   (5)\n|\t\t\tCt = C1||C2\t    (6)\n');

fprintf('|\n|\tThen, solving the system above by replacing (5) in (6),\n');
fprintf(fid,'|\n|\tThen, solving the system above by replacing (5) in (6),\n');

fprintf('|\n|\t\t\t ------------------ \n|\t\t\t|  L  :%7.2f %cH  |\n|\t\t\t|  C_1:%7.2f %cF  |\n|\t\t\t|  C_2:%7.2f %cF  |\n|\t\t\t ------------------ \n',L_adap,L_u,C_1,C_1_u,C_2,C_2_u);
fprintf(fid,'|\n|\t\t\t ------------------ \n|\t\t\t|  L  :%7.2f %cH  |\n|\t\t\t|  C_1:%7.2f %cF  |\n|\t\t\t|  C_2:%7.2f %cF  |\n|\t\t\t ------------------ \n',L_adap,L_u,C_1,C_1_u,C_2,C_2_u);

fprintf('|\n|\tFinally, we calculate the insertion loss of the network,\n');
fprintf(fid,'|\n|\tFinally, we calculate the insertion loss of the network,\n');

%% Calculate the network's insertion loss (IL)
%%% Get the IL equation
%%%% Impedances and other parameters definition
syms i1 i2 i3 Vs Vb Zb VL Zc1 Zc2 ZL RL RS i1_eq i2_eq i3_eq
Z1 = -1j/(w*C1);        % C1 capacity impedance  (Ω)
Z2 = -1j/(w*C2);        % C1 capacity impedancea (Ω)
Z3 = 1j*w*L;            % L  inductor impedance  (Ω)
%%%% Equivalent impedance Zb = ((RL||Z3)+Z2)||Z1
Zb    = 1./ ((1./((1./((1/Rl) + (1/Z3))) + Z2)) + (1/Z1));    % Numeric value (Ω)
Zb_eq = 1./ ((1./((1./((1/RL) + (1/ZL))) + Zc2)) + (1/Zc1)); % Symbolic value

Vb    = Vs*(Zb/(Rs+Zb));       % Numeric  value (V)
Vb_eq = Vs*(Zb_eq/(RS+Zb_eq)); % Symbolic value
%%%% Currents
i1    = (Vs-Vb)/Rs;    % Numeric  value (A)
i1_eq = (Vs-Vb_eq)/RS; % Symbolic value

i2    = Vb/Z1;     % Numeric  value (A)
i2_eq = Vb_eq/Zc1; % Symbolic value

i3    = i1-i2;       % Numeric  value (A)
i3_eq = i1_eq-i2_eq; % Symbolic value

%%%% VL with the resonant circuit
VL_r    = Vb-(i3*Z2);        % Numeric  value (V)
VL_r_eq = Vb_eq-(i3_eq*Zc2); % Symbolic value

%%%% VL without the resonant circuit
VL_nr    = Vs*Rl/(Rl+Rs); % Numeric  value (V)
VL_nr_eq = Vs*RL/(RL+RS); % Symbolic value

%%%% Get the IL equation
IL     = simplify(VL_r/VL_nr);       % Numeric   value (linear)
IL_eqn = simplify(VL_r_eq/VL_nr_eq); % Symbolic value

%%%% Show the IL equation
fprintf('|\n|\t\tIL = 20log10(VL_resonant/VL_noResonant)\t(7)');
fprintf(fid,'|\n|\t\tIL = 20log10(VL_resonant/VL_noResonant)\t(7)');

fprintf('\n|\n|\t\tSolving (7)...');
fprintf(fid,'\n|\n|\t\tSolving (7)...');

fprintf('\n|\n|\t\tIL = 20log10(%s)\n',IL_eqn);
fprintf(fid,'\n|\n|\t\tIL = 20log10(%s)\n',IL_eqn);

IL_dB = 20*log10(IL);
fprintf('|\n|\n|\t\t\t ------------------- \n|\t\t\t|  IL(dB) :%7.2f  |\n|\t\t\t ------------------- \n',IL_dB);
fprintf(fid,'|\n|\n|\t\t\t ------------------- \n|\t\t\t|  IL(dB) :%7.2f  |\n|\t\t\t ------------------- \n',IL_dB);

fprintf('*******************************************************************************************\n\n\n');
fprintf(fid,'*******************************************************************************************\n\n\n');


fclose(fid); %Close file

function [num_adap,unit] = units(orig_num)

if orig_num>=1e9
    num_adap      = orig_num*1e-9;
    unit = 'G';
elseif orig_num>=1e6
    num_adap      = orig_num*1e-6;
    unit = 'M';
elseif orig_num>=1e3
    num_adap      = orig_num*1e-3;
    unit = 'k';
elseif orig_num>=1e1
    num_adap      = orig_num;
    unit = '';
elseif orig_num<=1e-9
    num_adap      = orig_num*1e12;
    unit = 'p';
elseif orig_num<=1e-6
    num_adap      = orig_num*1e9;
    unit = 'n';
elseif orig_num<=1e-3
    num_adap      = orig_num*1e6;
    unit = 'u';
elseif orig_num<=1e-1
    num_adap      = orig_num*1e3;
    unit = 'm';
end

