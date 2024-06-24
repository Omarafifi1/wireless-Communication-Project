City_area = 100;
SIR_mindB = 14;
SIR_min = 10^(SIR_mindB/10);
user_dens = 100:2000; 
S = 340;                     %number of channels/cluster
f = 900;                     %frequency in MHz
h_bs = 20;                   %effective height of BS
h_ms = 1.5;                  %effective height of MS
Ms_sens = -95;               %MS sensitivity in dBm
n = 4;                       %path loss exponent
Au = 0.025;                  %traffic intensity per user in Erlang

GOS = 0.02;

for Secortiziton = [6 2 1]
    
if Secortiziton == 6
    n_sectors = 1;
elseif Secortiziton == 2
        n_sectors = 3;
elseif Secortiziton == 1
        n_sectors = 6;
end
%====================================================================================

N_min = ((((Secortiziton)*SIR_min)^(1/n) + 1)^2) / 3;   

possible_N = [];
for i = 0:10                          % making an array with i^2 + ik + k^2 values
    for k = 1:10
        all_N = i^2 + i*k + k^2;
        possible_N(end +1) = all_N;
    end
end


possible_N = sort(possible_N);
possible_N = unique(possible_N(:).');
for i = 1:length(possible_N)
    if N_min > possible_N(1,i) && N_min < possible_N(1,i+1)
        N = possible_N(1,i+1);                %Cluster size
        break
    end
end
%==========================================================================================

k = floor(S/(N*n_sectors));                                %channels/sector

syms A
eqn = ((A^k) / factorial(k)) == GOS * (sum(A.^(0:k) ./ factorial(0:k)));
A_sol = solve(eqn,A);
A_total = double(A_sol);

if length(A_total) > 1
    for i= 1:length(A_total)
        if isreal(A_total(i))
            A_sector = (A_total(i));                   %traffic/sector
        end
    end
end

A_cell = A_sector*n_sectors;                                %traffic/cell
%=========================================================================================
% calculating the total number of cells:

all_users = City_area .* user_dens;
A_network = Au .* all_users;
total_cells = ceil((A_network ./ A_cell));                     %total number of cells
%=========================================================================================
% calculating the cell radius of cell:

R = sqrt((2*(City_area./total_cells))/(3*sqrt(3)));

%=========================================================================================
%plotting:

if Secortiziton == 6
figure(1)
plot(user_dens,total_cells)
title('The number of cells versus user density , SIR = 14 dB','FontSize',14,'FontWeight','bold')
xlabel('User density','FontSize',12,'FontWeight','bold','Color','r')
ylabel('Number of cells','FontSize',12,'FontWeight','bold','Color','b')

figure(2)
plot(user_dens,R)
grid on
title('Cell Radius versus user density , SIR = 14 dB','FontSize',12,'FontWeight','bold')
xlabel('User density','FontSize',12,'FontWeight','bold','Color','r')
ylabel('Cell Radius','FontSize',12,'FontWeight','bold','Color','b')

elseif Secortiziton == 2
     figure(1)
     hold on
     plot(user_dens,total_cells)
     grid on
     
     figure(2)
     hold on
     plot(user_dens,R)
     grid on

elseif Secortiziton == 1
    figure(1)
    hold on
    plot(user_dens,total_cells)
    grid on

    figure(2)
    hold on
    plot(user_dens,R)
    grid on

end
end

legend('Omni directional','120 sectorization','60 sectorization')