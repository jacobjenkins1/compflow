%Compressible Flow Capstone Project
%Jacob E Jenkins
%Problem 1
clear; close; clc;

%Inputs
gamma = 1.4; % Ratio of Specific Heats
M1 = [2]; % Incoming Mach number
beta_1 = -40;

%Determine if the given beta is a left or right running wave
if beta_1 > 0
    beta_1_type = 'Left Running Wave';
else
    beta_1_type = 'Right Running Wave';
end

beta_1 = abs(beta_1); %Put beta back into a positive value for ease of theta calculations

n_M = length(M1);
figure()
theta = cell(1,n_M);

max_theta = zeros(1,n_M);
max_beta = zeros(1,n_M);
for i = 1:n_M
    M = M1(i);
    beta = asin(1/M):0.001:pi/2;
    theta{i} = atan( 2 .* cot(beta) .* (M^2 .* sin(beta).^2 - 1) ./(M^2 .* (gamma + cos(2.*beta)) + 2) );
    % Find max theta & beta
    [max_theta(i), max_idx] = max( theta{i} );
    max_beta(i) = beta(max_idx);
    plot(theta{i}.*180/pi, beta.*180/pi,'Color', 0.7*(1 - i/ n_M).*[1 1 1],'LineWidth', 1.5); hold on
    xlim([0,50]);
    ylim([0,90]);
end


theta_deg = theta{1,1}*180/pi;
beta_deg = beta*180/pi;
theta_beta = [theta_deg' beta_deg'];
theta_1 = theta_beta(theta_beta(:,2) > (beta_1 - 0.1) & theta_beta(:,2) < (beta_1 + 0.1),:);
theta_1 = mean(theta_1(:,1));




plot(max_theta.*180/pi, max_beta.*180/pi, '--k');
xlabel('Deflection Angle, $\theta$ degrees', 'Interpreter', 'Latex', 'fontsize',18);
ylabel('Shock wave angle, $\beta$ degrees','Interpreter', 'Latex', 'fontsize',18);
set(gca,'TickLabelInterpreter','latex')
grid on
scatter(theta_1,beta_1)
ax = gca;
ax.FontSize = 16;
hold off

M1_n = M1*sind(beta_1);
M2_n = sqrt((1+((gamma-1)/2)*M1_n^2)/((gamma*M1_n^2)-(gamma-1)/2));
M2 = M2_n/sind(beta_1-theta_1);
p2_p1 = 1+(((2*gamma)/(gamma+1))*(M1_n^2-1));
 
%%
M2 = [M2]; % Incoming Mach number
theta_2 = theta_1

n_M_2 = length(M2);
figure()
theta = cell(1,n_M_2);

max_theta_2 = zeros(1,n_M_2);
max_beta_2 = zeros(1,n_M_2);
for i = 1:n_M_2
    M = M2(i);
    beta = asin(1/M):0.001:pi/2;
    theta{i} = atan( 2 .* cot(beta) .* (M^2 .* sin(beta).^2 - 1) ./(M^2 .* (gamma + cos(2.*beta)) + 2) );
    % Find max theta & beta
    [max_theta(i), max_idx] = max( theta{i} );
    max_beta(i) = beta(max_idx);
    plot(theta{i}.*180/pi, beta.*180/pi,'Color', 0.7*(1 - i/ n_M).*[1 1 1],'LineWidth', 1.5); hold on
    xlim([0,50]);
    ylim([0,90]);
end

theta_deg = theta{1,1}*180/pi; % Convert turn angle from rad to deg
beta_deg = beta*180/pi; % Convert wave angle from rad to deg

theta_beta = [theta_deg' beta_deg']; % Array of theta with corresponding beta vals

% Cut out strong shock solutions
[val, idx] = max(theta_beta(:,1)); 
theta_beta = theta_beta(1:idx, :);

% Find the theta value which is closest to the input theta_2
theta_2 = theta_beta(theta_beta(:,1) > (theta_2 - 0.1) & theta_beta(:,2) < (theta_2 + 0.1),:);
theta_2 = mean(theta_1(:,1));
beta_2 = mean(theta_beta(:,2));

% Calculate wave angle with respect to the wall
phi_2 = beta_2 - theta_2;


M2_n = M2*sind(beta_2);
M3_n = sqrt((1+((gamma-1)/2)*M2_n^2)/((gamma*M2_n^2)-(gamma-1)/2));
M3 = M3_n/sind(beta_2-theta_2);
p3_p2 = 1+(((2*gamma)/(gamma+1))*(M2_n^2-1));
p3_p1 = p3_p2*p2_p1;




sprintf("Problem 1 Answers: \na) M_2 = %d \nb) phi = %f degrees \nc) M_3 = %f \nd) P3/P1 = %f ", M2, phi_2, M3, p3_p1)
