fc=900;  %%mhz
s=340;
h_MS=1.5;
h_bs=20;
G_MS=-95; %%dbm
n=4;
Au=0.025;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GOS = input("Enter the GOS : ");
city_area = input("Enter the city area : ");
user_density = input("Enter the user density : ");
min_SIR_db = input("Enter the min_SIR_db : ");
sectorization_method = input("Enter the sectorization_method press 6 for omnidirectional or 2 for 120 degrees or 1 for 60 degrees: ");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%          calaculate cluster_size               %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
min_SIR=power(10,(min_SIR_db/10));
y=power((min_SIR*sectorization_method),(1/n));
reuse_ratio = y+1;
z=power(reuse_ratio,2);
Nmin=z/3;

N_values=[];
for i=0:10
  for k=1:10
    N_values(i+1,k)=i^2+k^2+i*k;
  end
end
N_values=unique(N_values);

for j=1:(length(N_values))
  if(Nmin <= N_values(j))
  cluster_size=N_values(j);
  break
  end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%        calaculate number of sectors                %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sectorization_method ==1;
  number_of_sectors=6;
elseif sectorization_method==2;
  number_of_sectors=3;
elseif sectorization_method==6;
  number_of_sectors=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%        calaculate number of channels per sector    %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ch_per_sector=floor(s/(cluster_size*number_of_sectors));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%         calaculate traffic intensity per sector    %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms A
eqn = ((A^ch_per_sector) / factorial(ch_per_sector)) == GOS * (sum(A.^(0:ch_per_sector) ./ factorial(0:ch_per_sector)));
A_sol = solve(eqn,A);
A_total = double(A_sol);

if length(A_total) > 1
    for i= 1:length(A_total)
        if isreal(A_total(i))
            A_sector = round(A_total(i));                   %traffic/sector
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%         calaculate traffic intensity per cell      %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

traffic_intensity_per_cell=number_of_sectors*A_sector;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%         calaculate total_no_of_cells               %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_traffic = city_area * user_density * Au;
total_no_of_cells= ceil(total_traffic / traffic_intensity_per_cell);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%          calaculate cell_radius                         %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
area_of_hexagonal_cell=city_area/total_no_of_cells;
cell_radius=sqrt((2*area_of_hexagonal_cell)/(3*sqrt(3)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%          calaculate transmitted power          %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ch=0.8+(1.1*log10(fc)-0.7)*h_MS-1.56*log10(fc);
FSL=69.55+(26.16*log10(fc))-(13.82*log10(h_bs))-ch+((44.9-(6.55*log10(h_bs)))*log10(cell_radius));
P_TX=FSL+G_MS;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%          calaculate recieved power             %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting MS received power vs. distance from BS
R_range = [0:0.01:cell_radius];
P_rx = P_TX - (69.55 + 26.16*log10(fc) -13.82*log10(h_bs) - ch +(44.9 - 6.55*log10(h_bs))*log10(R_range));
plot(R_range,P_rx)
grid on
title('The MS received power in dBm versus the receiver distance from the BS')
xlabel('Distance(km)')
ylabel('Power-Rx(dBm)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%                % Displaying results            %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Cluster Size: %d\n', cluster_size);
fprintf('Number of Cells: %d\n', total_no_of_cells);
fprintf('Cell Radius: %.2f km\n', cell_radius);
fprintf('Traffic Intensity per Cell: %.2f Erlang\n', traffic_intensity_per_cell);
fprintf('Traffic Intensity per Sector: %.2f Erlang\n', A_sector);
fprintf('Base Station Transmitted Power: %.2f dBm\n', P_TX);
