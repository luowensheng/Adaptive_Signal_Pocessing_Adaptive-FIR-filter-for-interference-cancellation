

clc;
clear all
close all
tic
N = 50000; % number of samples
g = 0.5; % leakage factor
FL = 12; %filter length
m=100; % constant to make vizualizing plots clearer
n1=[0:1:N-1]; % initialize time
i=((1/3)*cos(2*pi*n1/3)+(1/4)*cos(2*pi*n1/4));
h=[0.5484 -0.4888 -0.4356 -0.3882 -0.3461]';
a= [0.01 0.0013 0.005];  %different alpha values
c = 0.001;
error = zeros(3,N);

W=zeros(3,3); % save the averaged |e(n)-s(n)|^2 for the last 5000 iterations
Z=zeros(3,3); % save the averaged |e(n)-s(n)|^2 after 100000 iterations



for K=1:length(a) % this loop iterates to test different values of alpha
    alpha= a(K);
    S=zeros(3,N); % save values for |e(n)-s(n)|^2
    error = zeros(3,N);
for trials=1:200

s=[ones(1,length(n1)/2)  (-1)*ones(1,length(n1)/2) ];
s=s(randperm(length(s)));  % s(n) vector


x= conv(i,h); 
x = x(1:N); % the reference input x(n)
d = s + i; %desired signal

%set/reset all fn to 0 for the update equations
f_LMS = zeros(1,FL)';
f_NLMS = zeros(1,FL)';
f_MLMS = zeros(1,FL)';
f_MLMS_pre = zeros(1,FL)';


%LMS Algorithm
for n=1:N
   if n < FL
       x_LMS = [x(1,n:-1:1) zeros(1,FL-n)].';
   else
       x_LMS = x(1,n:-1:n-FL+1).';
   end
   y_LMS = f_LMS.'*x_LMS;
   e_n(1,n) = d(n) - y_LMS;
   f_LMS = f_LMS + alpha*e_n(1,n)*x_LMS;
end
%NLMS  Algorithm
for n=1:N
   if n < FL
       x_NLMS = [x(1,n:-1:1) zeros(1,FL-n)].';
   else
       x_NLMS = x(1,n:-1:n-FL+1).';
   end
   y_NLMS = f_NLMS.'*x_NLMS;
   e_n(2,n) = d(n) - y_NLMS;
   f_NLMS = f_NLMS + (alpha/(c+x_NLMS.'*x_NLMS))*e_n(2,n)*x_NLMS;
end


%MLMS Algorithm
for n=1:N
   if n < FL
       x_MLMS = [x(1,n:-1:1) zeros(1,FL-n)].';
   else
       x_MLMS = x(1,n:-1:n-FL+1).';
   end
   f_MLMS_pre = f_MLMS;
   y_MLMS = f_MLMS.'*x_MLMS;
   e_n(3,n) = d(n) - y_MLMS;
   f_MLMS = f_MLMS + (1-g)*(f_MLMS-f_MLMS_pre)+alpha*g*(e_n(3,n))*x_MLMS;
end

error(1,:)=error(1,:)+(abs(e_n(1,:)).^2);   % save error for the LMS filter
error(2,:)=error(2,:)+(abs(e_n(2,:)).^2);   % save error for the NLMS filter
error(3,:)=error(3,:)+(abs(e_n(3,:)).^2);   % save error for the MLMS filter


S(1,:)=abs(e_n(1,:)-s).^2; %save |e(n)-s(n)|^2 for the LMS filter
S(2,:)=abs(e_n(2,:)-s).^2; %save |e(n)-s(n)|^2 for the NLMS filter
S(3,:)=abs(e_n(3,:)-s).^2; %save |e(n)-s(n)|^2 for the MLMS filter

end
%find average error
error_lms=error(1,:)/trials;
error_nlms=error(2,:)/trials;
error_mlms=error(3,:)/trials;

%plot the graphs
figure(K)
subplot(3,1,1),plot(error_lms(1:m:end))
ylabel('LMS error')
xlabel('Iterations')
title(sprintf(' Plot of of the averaged |e(n)|^2 for the LMS error alpha = %s',alpha))


subplot(3,1,2), plot(error_nlms(1:m:end))
ylabel('NLMS error')
xlabel('Iterations')
title(sprintf(' Plot of of the averaged |e(n)|^2 for the NLMS error alpha = %s',alpha))


subplot(3,1,3),plot(error_mlms(1:m:end))
ylabel('Momentum LMS error')
xlabel('Iterations')
title(sprintf(' Plot of of the averaged |e(n)|^2 for the Momentum LMS error alpha = %s',alpha))




% Find the averaged |e(n)-s(n)|^2 for the last 5000
W(1,K)=mean(S(1,end-5000:end));   % for LMS filter
W(2,K)=mean(S(2,end-5000:end));   % for NLMS filter
W(3,K)=mean(S(3,end-5000:end));   % for MLMS filter


% Find the averaged |e(n)-s(n)|^2 after 100000 iterations
Z(1,K)=mean(S(1,1:10000));      % for LMS filter
Z(2,K)=mean(S(2,1:10000));      % for NLMS filter
Z(3,K)=mean(S(3,1:10000));      % for MLMS filter





end
display('a= 0.01      0.0013     0.005')
display('   0.01      0.0013     0.005')
 display('The averaged |e(n)-s(n)|^2 for the last 50000 respectively:')
display(W)

display(' The averaged |e(n)-s(n)|^2 for the last 100000 iterations:')
display(Z)
toc


