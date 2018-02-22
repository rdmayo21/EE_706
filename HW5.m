%%
%Problem 1

% Scalar Kalman Filter Example
% Copyright 2018 M. C. Budge Jr.
clear;close all
kup=50;
kplt=[0:kup];ksamp=length(kplt)-1;

% Set up a Monte Carlo loop
% Monte Carlo on xa and xh
Monte = 1; % Number of Monte Carlo runs
% Initialize the Monte Carlo arrays
xamonte = zeros(Monte,kup+1);
xhmonte = xamonte;

for Q = 0:0.00666:.1
    for monte = 1:Monte

        %%%%%%%%%%generate the actual state and the measurement %%%%%%%%
        Ra=0.01; % actual measurement variance
        Qa=0; % Variance on actual system disturbance; setting to zero for hw
        xa=zeros(size(kplt)); % actual system state
        ya=xa;	% set the array size
        xa(1)=1; % The initial value is 1.0
        ya(1)=xa(1); % Set the initial measurement to
        for k=1:ksamp
            xa(k+1)=0.97*xa(k)+sqrt(Qa)*randn(1);
            ya(k+1)=xa(k+1)+randn(1)*sqrt(Ra);
        end
        xamonte(monte,:) = xa;
        %  **** DEFINE THE MODEL PARAMETERS *****
        F=0.9; % state transition matrix
        G=1;		% input distribution matrix
        H=1;		% output allocation matrix
%         Q=1; % disturbance variance
        R=0.01; % measurement variance

        % zero all of the Kalman filter parameters
        xp=zeros(size(xa)); % predicted state xh(k+1|k)
        xh=xp; % state estimate
        K=xp; % Kalman gain
        Pp=xp; % predicted covariance P(k+1|k)
        yp=xp; % predicted measurement yh(k+1|k)
        P=xp; % State Covariance
        xh(1)=0; % Set the initial value of the state estimate to 0
        P(1)=10; % Set the initial covariance to 10
        K(1)=1; % Kalman Gain = 1 since P(0) is so large
        %
        % ***** Implement the Kalman Filter *****
        for k=1:ksamp
            xp(k+1)=F*xh(k); % The predicted state
            Pp(k+1)=F*P(k)*F'+G*Q*G'; % The predicted covariance
            K(k+1)=H*Pp(k+1)/(H*H*Pp(k+1)+R); % The Kalman gain
            yp(k+1)=H*xp(k+1); % The predicted measurement
            xh(k+1)=xp(k+1)+K(k+1)*(ya(k+1)-yp(k+1)); % The smoothed state
            P(k+1)=(1-K(k+1))*Pp(k+1); % The smoothed covariance
        end
        xhmonte(monte,:) = xh;
    end
    %
    % Generate the plots
%     h1=figure(1);
    %First Plot the state and its estimate
%     sh1=subplot(221);
    figure
    plot(kplt,mean(xamonte,1),kplt,mean(xhmonte,1),'r','LineWidth',1.5)
    grid on
    title(['Q = ' num2str(Q)])
    xlabel('stage (k)')
    ylabel('value')
    legend('actual state','state estimate')
    % Next plot the Kalman gain
%     sh2=subplot(224);
%     plot(kplt,K,'LineWidth',1)
%     grid on
%     xlabel('stage (k)')
%     ylabel('K(k)')
%     zoom yon
%     % Next plot the Estimate Covariance
%     sh3=subplot(223);
%     plot(kplt,(P),'LineWidth',1.5)
%     grid on
%     ylabel('covariance (P)')
%     xlabel('stage (k)')
%     zoom yon
%     % Next plot the smoothed rms error (sqrt(P))
%     sh3=subplot(222);
%     plot(kplt,sqrt(P),'LineWidth',1.5)
%     grid on
%     ylabel('rms error (sqrt(P))')
%     xlabel('stage (k)')
%     zoom yon
end
% set(h1,'Position',[360 133 559 565])


%% 
%Problem 2

clear;close all
kup=50;
kplt=[0:kup];ksamp=length(kplt)-1;

% Set up a Monte Carlo loop
% Monte Carlo on xa and xh
Monte = 1; % Number of Monte Carlo runs
% Initialize the Monte Carlo arrays
xamonte = zeros(Monte,kup+1);
xhmonte = xamonte;

for Q = 0:0.1:1
    for monte = 1:Monte

        %%%%%%%%%%generate the actual state and the measurement %%%%%%%%
        Ra=0.01; % actual measurement variance
        Qa=0; % Variance on actual system disturbance; setting to zero for hw
        xa=zeros(size(kplt)); % actual system state
        ya=xa;	% set the array size
        xa(1)=1; % The initial value is 1.0
        ya(1)=xa(1); % Set the initial measurement to
        for k=1:ksamp
            xa(k+1)=0.97*xa(k)+sqrt(Qa)*randn(1);
            ya(k+1)=xa(k+1)+randn(1)*sqrt(Ra);
        end
        xamonte(monte,:) = xa;
        %  **** DEFINE THE MODEL PARAMETERS *****
        F=0.97; % state transition matrix
        G=1;		% input distribution matrix
        H=0.8;		% output allocation matrix
%         Q=1; % disturbance variance
        R=0.01; % measurement variance

        % zero all of the Kalman filter parameters
        xp=zeros(size(xa)); % predicted state xh(k+1|k)
        xh=xp; % state estimate
        K=xp; % Kalman gain
        Pp=xp; % predicted covariance P(k+1|k)
        yp=xp; % predicted measurement yh(k+1|k)
        P=xp; % State Covariance
        xh(1)=0; % Set the initial value of the state estimate to 0
        P(1)=10; % Set the initial covariance to 10
        K(1)=1; % Kalman Gain = 1 since P(0) is so large
        %
        % ***** Implement the Kalman Filter *****
        for k=1:ksamp
            xp(k+1)=F*xh(k); % The predicted state
            Pp(k+1)=F*P(k)*F'+G*Q*G'; % The predicted covariance
            K(k+1)=H*Pp(k+1)/(H*H*Pp(k+1)+R); % The Kalman gain
            yp(k+1)=H*xp(k+1); % The predicted measurement
            xh(k+1)=xp(k+1)+K(k+1)*(ya(k+1)-yp(k+1)); % The smoothed state
            P(k+1)=(1-K(k+1))*Pp(k+1); % The smoothed covariance
        end
        xhmonte(monte,:) = xh;
    end
    %
    % Generate the plots
%     h1=figure(1);
    %First Plot the state and its estimate
%     sh1=subplot(221);
    figure
    plot(kplt,mean(xamonte,1),kplt,mean(xhmonte,1),'r','LineWidth',1.5)
    grid on
    title(['Q = ' num2str(Q)])
    xlabel('stage (k)')
    ylabel('value')
    legend('actual state','state estimate')
    % Next plot the Kalman gain
%     sh2=subplot(224);
%     plot(kplt,K,'LineWidth',1)
%     grid on
%     xlabel('stage (k)')
%     ylabel('K(k)')
%     zoom yon
%     % Next plot the Estimate Covariance
%     sh3=subplot(223);
%     plot(kplt,(P),'LineWidth',1.5)
%     grid on
%     ylabel('covariance (P)')
%     xlabel('stage (k)')
%     zoom yon
%     % Next plot the smoothed rms error (sqrt(P))
%     sh3=subplot(222);
%     plot(kplt,sqrt(P),'LineWidth',1.5)
%     grid on
%     ylabel('rms error (sqrt(P))')
%     xlabel('stage (k)')
%     zoom yon
end
% set(h1,'Position',[360 133 559 565])
