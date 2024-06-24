City_area = 100;     %Km square 
%=======================================================
x=1;
for sectroization = [6 2 1] %%%%%%%%%%%%%%%%%% for loop to plot the 3 cases
SIR_range = 1:0.01:30;

% making an array of cluster size for different SIR:
N_min = [];
for SIR_min = 10.^(SIR_range/10) %%%%%%%%  db to ratio
N_i = ((((sectroization)*SIR_min)^(1/4) + 1)^2) / 3;
N_min(end+1) = N_i; %%% appends the calculated value N_i to the end of the array N_min
end

possible_N = [];  %%%calculate the possible values of N for different i,k
for i = 0:10                          % making an array with i^2 + ik + k^2 values
    for k = 1:10
        all_N = i^2 + i*k + k^2;
        possible_N(end +1) = all_N;
    end
end
possible_N = sort(possible_N);  %%%%%%%% sorts the elements of possible_N in ascending order.
possible_N = unique(possible_N(:).');%%%%%%%% removes duplicate elements from

j = 1;
for outer = 1:length(N_min) 
for i = 1:length(possible_N)
    if N_min(j) > possible_N(1,i) && N_min(j) < possible_N(1,i+1) %%%%%%%% the cluster size final values
        N(j) = possible_N(1,i+1);                %Cluster size
        j = j+1;
        break
    end
    
end
end

if sectroization == 6       %%%%%%%case 1
    figure(1)
    grid on
    hold on
    plot(SIR_range,N)
    title('The cluster size versus SIR with range  1dB to 30dB','FontSize',14,'FontWeight','bold','color','b')
    xlabel('SIR-min (dB)')
    ylabel('Cluster size (N)')
    x=x+1;
elseif sectroization == 2       %%%%%%%%%%case2
    hold on
    plot(SIR_range,N)
    
elseif sectroization == 1       %%%%%%%%%%case3
    hold on
    plot(SIR_range,N)
    
end
end
%========================================================================
legend('Omni directional','120 sectroization','60 sectroization')