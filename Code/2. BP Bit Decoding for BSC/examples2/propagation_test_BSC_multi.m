%test and plot
%J=10 and Q=200 approximately runs for 2 mins 


pcm=[0 0 1 0 0 1 1 1 0 0 0 0 ; 1 1 0 0 1 0 0 0 0 0 0 1; ...
    0 0 0 1 0 0 0 0 1 1 1 0 ; 0 1 0 0 0 1 1 0 0 1 0 0; ...
    1 0 1 0 0 0 0 1 0 0 1 0 ; 0 0 0 1 1 0 0 0 1 0 0 1; ...
    1 0 0 1 1 0 1 0 0 0 0 0 ; 0 0 0 0 0 1 0 1 0 0 1 1; ...
    0 1 1 0 0 0 0 0 1 1 0 0];
J=10;% how many error values are tested
Q=200;% number of codewords sent each test
upper=0.2;

es=linspace(upper/J,upper-upper/J,J-1); %error values

BERS=zeros(J-1,3); %BER for every algorithm
times=zeros(J-1,3); %time for every algorithm

for i=1:J-1
    [BERS(i,:),times(i,:)]= propagation_test_BSC_f(pcm,es(i),Q);
end

figure
hold on
plot(es,BERS(:,1))
plot(es,BERS(:,2))
plot(es,BERS(:,3))
plot (es,es)
set(gca, 'YScale', 'log')
legend('bitDecode','decode','Map Decoding','e')
title('Bit error rate for various algorithms')
xlabel('Channel e parameter')
ylabel('BER')
figure
hold on
plot(es,times(:,1))
plot(es,times(:,2))
plot(es,times(:,3))

