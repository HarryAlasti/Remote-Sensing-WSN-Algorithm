clear all; clc;
close all;

L = 10000;                                                    % # of randomly distributed sensors in the field
n_std = 0.1;                                                  % random wite noise std in sensor observation
Nlevels0 = 2;                                                % Initial # of contour levels in spatial modeling
lambda = 0.005;                                           % tolerance bound of convergence factor
Rough = 0.25;                                                % ratio of random Gaussian 2D dist in remotr sensing (rs) spatial modeling
jmp = 1;                                                        % the contour line increment in modeling spatial in initiation phase
Npnt = 10;                                                     % a factor of  updates in temporal phase
grad_type = 2;                                               % grad_type = 1: for type-1 adapted gradient (agressive)
                                                                      % grad_type = 2: for type-2 adapted gradient 
                                                                      % grad_type = 0: no gradient
                                                                      
LVL_Type = 1;                                               % LVL_Type = 0 mean equaly-spaced levels, 
                                                                       % LVL_Type = 1 means Lloyd-Max
kap_Flag = 1;                                                 % Kappa factor = 1 :unequal level-increment
                                                                       
MA =10;                                                         % moving average window size
dyn_end = 300;                                               % total temporal instances
Step = 0.001;                                                  % amount of Gaussian center shifts in temporal (for x and y)
% =================================================================
%     Spatial Modelling using Wireless Sensing
% =================================================================
% memory allocation for the used blocks
Sensors = zeros(L,4);

%% Loading the saved spatial distribution's data
name = sprintf('NewWave_%d',1);
load(name);

V = 0:1:100;
[XI,YI] = meshgrid(V);
Zrs = zeros(size(XI));
for k = 1: round(N * Rough)
    Zrs = Zrs +Peak * (2 * (P(k,3) - 0.5))/ sqrt(2*pi) / Sigma * exp(-((XI-P(k,1)).^2 + (YI-P(k,2)).^2)/2/Sigma^2);
end

% Snapshot normalization
Zrs = Peak / 2 * Zrs/(max(max(Zrs)) - min(min(Zrs))) + 50;
Data_rs = Zrs;
% Plot of the actual, remote sensing spatial signal, and the difference
%{
figure(18); mesh(XI,YI,Data); hold on;
Vec = 35:5:70;
figure(18); contour3(XI,YI,Data,Vec);

figure(19); mesh(XI,YI,Data_rs); hold on;
figure(19); contour3(XI,YI,Data_rs,Vec);

Error = abs(Data-Data_rs);
figure(20); mesh(XI,YI,Error); hold on;
%}
%% Defining L uniformly distributed wireless sensors over the area 
Sensors(:,1:2) = Peak * rand(L,2);                 % Defines a random distribtion of Ns sensors
X = Sensors(:,1);
Y = Sensors(:,2);
Z = zeros(L,1);
for k = 1: round(N * Rough)
    Z = Z + 2*Peak * (P(k,3) - 0.5) / sqrt(2*pi) / Sigma * exp(-((X-P(k,1)).^2 + (Y-P(k,2)).^2)/2/Sigma^2);
end

% Snapshot normalization
Z = Peak / 2 * Z / (max(max(Z)) - min(min(Z)))+50;
Sensors(:,3) = Z;        
% Plot & verification of the sensor mappings
%{
figure(1); plot3(Sensors(:,1),Sensors(:,2),Sensors(:,3),'.'); 
ZI = griddata(Sensors(:,1),Sensors(:,2),Sensors(:,3),XI,YI,'v4');          % linear   nearest      cubic     v4(biharmonic interpolation)
figure(2); mesh(XI,YI,abs(ZI-Data));
pause();
%}
%%
Nlevels = Nlevels0;
L0 = min(min(Zrs));
Lend = max(max(Zrs));

Error = [];
Error0 = [];
LEVELS = [Nlevels0];
Delta = [];
gamma = [0 , 0];

rpt = 0;
Flag = 0;
Nsense = 0;
Data_selkt = []; 
while (Flag == 0)
    rpt = rpt + 1;
   % Data_selkt = []; 

    if (LVL_Type == 0)
        Levels = linspace(L0,Lend,Nlevels+2);
        Levels = Levels(2:end-1);
    elseif (LVL_Type == 1)        
        [p,xx] = hist(Z(:),500);
        Levels = Lloyd_Max_2(p, xx, Nlevels);
    end
    Diff_lvl = abs(Levels(1:end-1) - Levels(2:end));
    Delta = [Delta , 0.2 * min(Diff_lvl)];

    for cnt = 1 : length(Levels)
        lvl = Levels(cnt);
        Npick = find(abs(Sensors(:,3) - lvl) < Delta(rpt));
        Data_selkt = [Data_selkt; Sensors(Npick,:)];
        Nsense = Nsense + length(Npick);
    end    
     
    Znew = griddata(Data_selkt(:,1), Data_selkt(:,2), Data_selkt(:,3), XI, YI,'v4');  
    error = mean(abs(Znew(:) - Zrs(:)));
    Error = [Error , error];    
    error = mean(abs(Znew(:) - Data(:)));
    Error0 = [Error0 , error];  

     if (rpt > 2)
         xx = 3* Error(rpt) / (Error(rpt) + Error(rpt-1) + Error(rpt-2)) ;
         gamma = [gamma , xx];
         if ((xx > 1 - lambda) && (xx < 1+lambda))
             Flag = 1;
         end
     end

    Nlevels = Nlevels + 1;
    if (Flag == 0)
        LEVELS = [LEVELS , Nlevels];
    end   
end
%{
figure(4); plot(LEVELS , 20*log10(Error),'b'); xlabel('# of Levels'); ylabel('MAE (dB)'); grid on; hold on;
figure(5); plot(LEVELS , Delta,'b');  xlabel('# of Levels'); ylabel('\Delta');grid on; hold on;
figure(6); plot(gamma); xlabel('# of iterations (n)'); ylabel('\gamma_n'); grid on;
%}
Delta0 = Delta(end);
M = LEVELS(end);
SNS = Sensors(:,1:2);
save('Initial_phase','Levels','M','Delta0','L0','Lend','SNS', 'L','lambda');

% =================================================================
%     Spatial Modelling using Wireless Sensing
% =================================================================
ZI = zeros(size(XI));
for k = 1: N
    ZI = ZI + 2*Peak * (P(k,3) - 0.5) / sqrt(2*pi) / Sigma * exp(-((XI-P(k,1)).^2 + (YI-P(k,2)).^2)/2/Sigma^2);
end
% Genmeration of small scale Gaussians
for k = 1 : Ns
    ZI = ZI + 5 * (Ps(k,3) - 0.5) / sqrt(2*pi) / Ss * exp(-((XI-Ps(k,1)).^2 + (YI-Ps(k,2)).^2)/2/Ss^2);
end
% Snapshot normalization
ZI = Peak / 2 * ZI / (max(max(ZI)) - min(min(ZI))) + 50;

% Plot of the actual, remote sensing spatial signal, and the difference
%{
figure(18); mesh(XI,YI,Data); hold on;
Vec = 35:5:70;
figure(18); contour3(XI,YI,Data,Vec);

figure(19); mesh(XI,YI,Data); hold on;
figure(19); contour3(XI,YI,Data,Vec);

Error = abs(Data-Data_rs);
figure(20); mesh(XI,YI,Error); hold on;
%}

%% Defining L uniformly distributed wireless sensors over the area 
Sensors(:,1:2) = SNS;                 % Defines a random distribtion of Ns sensors
X = Sensors(:,1);
Y = Sensors(:,2);
Z = zeros(L,1);
for k = 1: N
    Z = Z + 2*Peak * (P(k,3) - 0.5) / sqrt(2*pi) / Sigma * exp(-((X-P(k,1)).^2 + (Y-P(k,2)).^2)/2/Sigma^2);
end
% Genmeration of small scale Gaussians
for k = 1 : Ns
    Z = Z + 5 * (Ps(k,3) - 0.5) / sqrt(2*pi) / Ss * exp(-((X-Ps(k,1)).^2 + (Y-Ps(k,2)).^2)/2/Ss^2);
end
% Snapshot normalization
Z = Peak / 2 * Z / (max(max(Z)) - min(min(Z))) + 50;
Sensors(:,3) = Z;        
Sensors(:,4) = Z + n_std * randn(size(Z));          
% Plot & verification of the sensor mappings
%{
figure(1); plot3(Sensors(:,1),Sensors(:,2),Sensors(:,3),'.'); 
ZI = griddata(Sensors(:,1),Sensors(:,2),Sensors(:,3),XI,YI,'v4');          % linear   nearest      cubic     v4(biharmonic interpolation)
figure(2); mesh(XI,YI,abs(ZI-Data));
%}
%%
Nlevels = M;

Error2 = [];
Error2_lrn = [];
LEVELS2 = [M];
Delta2 = [Delta0];
SR = [];
gamma2 = [0 , 0];

rpt = 0;
Flag = 0;
Zold = Znew;
NSense = 0;
Data_selkt = []; 
while (Flag == 0)
    rpt = rpt + 1;
    %Data_selkt = []; 

    if (LVL_Type == 0)
        Levels = linspace(L0,Lend,Nlevels+2);
        Levels = Levels(2:end-1);
    elseif (LVL_Type == 1)        
        [p,xx] = hist(Z(:),500);
        Levels = Lloyd_Max_2(p, xx, Nlevels);
    end

    for cnt = 1 : length(Levels)
        lvl = Levels(cnt);
        Npick = find(abs(Sensors(:,4) - lvl) < Delta2(end));
        Data_selkt = [Data_selkt; Sensors(Npick,:)];
        NSense = NSense + length(Npick);
    end
    
    Znew = griddata(Data_selkt(:,1), Data_selkt(:,2), Data_selkt(:,4), XI, YI,'v4');  
    error = mean(abs(Znew(:) - ZI(:)));
    Error2 = [Error2 , error];    
    error =  mean(abs(Znew(:) - Zold(:)));
    Error2_lrn = [Error2_lrn , error];       
    
    %{
    if (rpt == 1)
        Diff_lvl = abs(Levels(1:end-1) - Levels(2:end));
        Delta2 = [Delta2 , 0.5 * min(Diff_lvl)];
    end
    %}
    SR = [SR , (max(Znew(:)) - min(Znew(:)))/(max(Zold(:)) - min(Zold(:)))];

    if (rpt > 1)        
       if (grad_type == 1)
            delta = Delta2(end) * SR(end) * (1 + (error - error_old)/(error + error_old));        
       elseif (grad_type == 2)                
            delta = Delta2(end) * (1 + (SR(end) / (1+SR(end)))^0.5 *(error - error_old)/(error + error_old));        
       elseif (grad_type == 0)
           delta = 0.2 * min(Levels(1:end-1) - Levels(2:end));
       end
       Delta2 = [Delta2 , delta];
       if (kap_Flag == 1)
                jmp = jmp_old + ceil(1 + 2*abs(Error2_lrn(end)-Error2_lrn(end-1))/(Error2_lrn(end)+Error2_lrn(end-1)));
       end
    end        

     if (rpt > 2)
         %xx = 3* Error2_lrn(rpt) / (Error2_lrn(rpt) + Error2_lrn(rpt-1) + Error2_lrn(rpt-2)) ;
          xx = 3* SR(rpt) / (SR(rpt) + SR(rpt-1) + SR(rpt-2)) ;
         gamma2 = [gamma2 , xx];
         if ((xx > 1 - lambda) && (xx < 1+lambda))
             Flag = 1;
         end
     end

    Nlevels = Nlevels + jmp;
    if (Flag == 0)
        LEVELS2 = [LEVELS2 , Nlevels];
    end  

    Zold = Znew;
    error_old = error;
    jmp_old = jmp;
end
Nlevels = Nlevels - jmp;

figure(14); plot(LEVELS , 20*log10(Error0),'k', LEVELS , 20*log10(Error),'b', LEVELS2 , 20*log10(Error2_lrn),'r', LEVELS2 , 20*log10(Error2),'m'); xlabel('# of Levels'); ylabel('MAE (dB)'); grid on; hold on;
%figure(14); plot(LEVELS , 20*log10(Error),'b', LEVELS2 , 20*log10(Error2_lrn),'r', LEVELS2 , 20*log10(Error2),'m'); xlabel('# of Levels'); ylabel('MAE (dB)'); grid on; hold on;
figure(15); plot(LEVELS , Delta,'b',LEVELS2 , Delta2,'r');  xlabel('# of Levels'); ylabel('\Delta');grid on; hold on;
figure(16); plot(gamma,'b'); xlabel('# of iterations (n)'); ylabel('\gamma_n'); grid on; hold on;
figure(16); plot(gamma2,'r'); xlabel('# of iterations (n)'); ylabel('\gamma_n'); grid on; hold on;
20*log10(Error2(end))
Involvement_percent = NSense / Nsense * 100
format short e;
[Nsense/L*100,  NSense/L*100,  NSense / Nsense * 100]
pause(1);
% =================================================================
%     Tracking 
% =================================================================

Error_opt = zeros(1, Npnt);
Error_unif = zeros(1, Npnt);
Error_learn_unif = zeros(1, Npnt);
Delta_unif_avg = zeros(1, Npnt);
Slevel = zeros(1, Npnt);
Sensors = zeros(L,7+MA+2);

Cost_opt = zeros(size(Slevel));
Cost_unif = zeros(size(Slevel));
Span_u = zeros(size(Slevel));
Slevel = zeros(2 , Npnt);

Cost = zeros(ceil(dyn_end/MA),1);
Error = zeros(ceil(dyn_end/MA),1);
COST = zeros(size(Cost));
ERROR = zeros(size(Error));

Delta_done = Delta2(end);

Sensors(:,1:2) = SNS;

x_s = Step;     % x direction shift increment
y_s = Step;     % y direction shift increment
reps = 0;
%Sensors(:,end) = Sensors(:,4);
Cost = 0 * Cost;
Error = 0 * Error;
m = 1;
Error_0 = Error_unif(end);
xx0 = 0;
yy0 = 0;
Delta = Delta2(end);
for dyn = 1:dyn_end
% Generation of the large scale Gaussians
    X = Sensors(:,1);
    Y = Sensors(:,2);
    Z = zeros(L,1);
    Data_new = zeros(size(XI));

   % if (abs(cos(2*pi*dyn*step)) > Thr)
        xx0 = xx0 + x_s * MA;
        yy0 = yy0 + y_s * MA;
        reps = reps + 1;
    %end

    for k = 1: N
        Z = Z + 2*Peak * (P(k,3) - 0.5) / sqrt(2*pi) / Sigma * exp(-((X-P(k,1)-xx0).^2 + (Y-P(k,2)-yy0).^2)/2/Sigma^2);
        Data_new = Data_new + 2*Peak * (P(k,3) - 0.5) / sqrt(2*pi) / Sigma * exp(-((XI-P(k,1)-xx0).^2 + (YI-P(k,2)-yy0).^2)/2/Sigma^2);
    end

% Generation of small scale Gaussians
    for k = 1 : Ns
        Z = Z + 5 * (Ps(k,3) - 0.5) / sqrt(2*pi) / Ss * exp(-((X-Ps(k,1)).^2 + (Y-Ps(k,2)).^2)/2/Ss^2);
        Data_new = Data_new + 5 * (Ps(k,3) - 0.5) / sqrt(2*pi) / Ss * exp(-((XI-Ps(k,1)).^2 + (YI-Ps(k,2)).^2)/2/Ss^2);
    end
    Z = Peak / 2 * Z / (max(Z) - min(Z))+50;
    Data_new = Peak / 2 * Data_new / (max(max(Data_new)) - min(min(Data_new)))+50;

    Sensors(:,4) = Z + n_std * randn(size(Z));     
    Sensors(:,5:7) = 0;

    if (dyn > MA)
        if (rem(dyn,MA) == 0)
            for clevel = 1 : length(Levels)
                index = find(abs(Sensors(:,4) - Levels(clevel)) <= Delta);  
                Sensors(index, 6) = 1;
            end
            Index = find(Sensors(:, 6) == 1);
            Cost(m) = length(Index); %+ length(diss);
            COST(m) = COST(m) + Cost(m);

            X1 = Sensors(Index, 1);
            Y1 = Sensors(Index, 2);
            Z1 = Sensors(Index, 4);
            Zx = griddata(X1,Y1,Z1,XI,YI,'v4');          % linear   nearest    cubic     v4(biharmonic interpolation)
            Error(m) = sum(sum(abs(Zx - Data_new))/(101*101));

            ERROR (m) = ERROR(m) + Error(m);
            m = m+1;

            L_0 = min(Zx(:));
            L_end = max(Zx(:));

     if (LVL_Type == 0)
        Levels = linspace(L_0,L_end, length(Levels)+2);
        Levels = Levels(2:end-1);
    elseif (LVL_Type == 1)        
        [p,xx] = hist(Zx(:),500);
        Levels = Lloyd_Max_2(p, xx, length(Levels));
    end

%{            
            Levels = linspace(L_0,L_end,Nlevels+2);
            Levels = Levels(2:end-1);
%}
        end
    end    
end

ERROR(1) = ERROR(2);
%figure(101); plot(20*log10(Error),'r*'); grid;title('performance'); ylabel('MAE in dB'); xlabel('update instant');
%figure(100); plot(Cost/L*100,'bs'); grid;title('cost'); ylabel('Percentage of sensors'); xlabel('update instant');

Len1 = (1:length(ERROR)) * MA;
figure(101); plot(Len1(1:end-1) , 20*log10(ERROR(1:end-1)),'r*'); grid;title('performance'); ylabel('MAE in dB'); xlabel('update instant');
figure(100); plot(Len1(1:end-1) , COST(1:end-1)/L*100,'bs'); grid;title('cost'); ylabel('Percentage of sensors'); xlabel('update instant');
%figure(320); mesh(XI,YI,Zx); grid on;
 

