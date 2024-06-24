City_area = 100;
SIR_mindB = 14;
SIR_min = 10^(SIR_mindB/10);   %%%%%db to ratio
user_density = 1400; 
S = 340;                     %number of channels/cluster
f = 900;                     %frequency in MHz
h_bs = 20;                   %effective height of BS
h_ms = 1.5;                  %effective height of MS
Ms_sens = -95;               %MS sensitivity in dBm
n = 4;                       %path loss exponent
Au = 0.025;                  %traffic intensity per user in Erlang
%=================================================================

for sectroization_method = [6 2 1]
   if sectroization_method == 6
    n_sectors = 1;
   elseif sectroization_method == 2
        n_sectors = 3;
   elseif sectroization_method == 1
        n_sectors = 6;
   end
   
GOS_range = 0.01:0.005:0.3;
N_min = ((((sectroization_method)*SIR_min)^(1/n) + 1)^2) / 3;   

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

k = floor(S/(N*n_sectors));                                %channels/sector
%====================================================================
q=1;
p=1;

% calculating traffic per cell for different GOS:
for GOS = 0.01:0.005:0.3  
syms A
eqn = ((A^k) / factorial(k)) == GOS * (sum(A.^(0:k) ./ factorial(0:k)));
A_sol = solve(eqn,A);
A_total = double(A_sol);

if length(A_total) > 1
    for i= 1:length(A_total)
        if isreal(A_total(i)) %%%%%to get the real values only
            A_sector = (A_total(i));                   %traffic/sector
        end
    end
end
A_cell = A_sector*n_sectors;
A_cell_array(p) = A_cell;                % vector of traffic per cell for every GOS
p = p + 1;
all_users = City_area * user_density;
A_network = Au * all_users;
total_cells(q) = ceil(A_network / A_cell);      
q=q+1;
end
%============================================================================
% plots:

if sectroization_method == 6
    figure(1)
    plot(GOS_range.*100,total_cells)
    title('The number of cells versus GOS (1% to 30%) at SIR = 19dB','FontSize',12,'FontWeight','bold','color','r')
    xlabel('GOS (%)')
    ylabel('Number of cells')

    figure(2)
    plot(GOS_range.*100,A_cell_array)
    title('The  traffic intensity per cell versus GOS (1%to 30%) at SIR=19dB','FontSize',12,'FontWeight','bold','color','r')
    xlabel('GOS (%)')
    ylabel('Traffic intensity per cell (A)')


elseif sectroization_method == 2
    figure(1) 
    hold on
    plot(GOS_range.*100,total_cells)

    figure(2)
    hold on
    plot(GOS_range.*100,A_cell_array)
        
elseif sectroization_method == 1
    figure(1) 
    hold on
    plot(GOS_range.*100,total_cells)
    grid on

    figure(2)
    hold on
    plot(GOS_range.*100,A_cell_array)
    grid on
end

end

legend('Omni directional','120 sectorization','60 sectorization')