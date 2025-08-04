%test and plot
%J=10 and Q=200 approximately runs for 2 mins 


%Designed pcm and gm
pcm =[1 1 1 1 0 0 1 0 0 0 0 0; 1 0 0 0 1 1 1 1 0 0 0 0;...
     0 1 1 1 0 0 0 1 1 0 0 0; 0 0 1 0 0 1 0 1 0 1 0 1; ...
     0 0 0 1 1 0 0 0 1 1 1 0; 1 0 0 0 0 0 0 0 0 0 1 1];
gm=[1 0 0 0 1 0 1 1 1 0 0 1; 0 0 1 1 1 1 0 0 0 1 1 1;...
    1 1 1 1 0 0 0 1 0 1 0 1; 0 1 0 0 0 1 1 0 1 0 1 1;...
    1 0 1 0 0 1 0 0 1 0 1 0; 0 0 0 1 0 0 1 1 0 0 1 1];

J=25;% how many error values are tested
Q=1000;% number of codewords sent each test
upper=0.3;

es=linspace(upper/J,upper-upper/J,J-1); %error values

BERS=zeros(J-1,3); %BERS for BP decoding vs MAP decoding

for i=1:J-1
    BERS(i,:)= propagation_test_BEC_f(pcm,gm,es(i),Q);
end

figure
hold on
plot(es,BERS(:,1),'linewidth',2)
plot(es,BERS(:,2),'linewidth',2)
plot(es,BERS(:,3),'linewidth',2)
plot (es,es)
set(gca, 'YScale', 'log')
legend('Bit decoding','Bit decoding wrong bit','Map Decoding','e')
title('Bit error rate comparison')
xlabel('Channel e parameter')
ylabel('BER')
grid on