%bit streams are character arrays here

%pcm = [ 1 1 1 1 0 0; 0 1 1 0 1 0; 1 0 1 0 0 1];
% pcm=[0 0 1 0 0 1 1 1 0 0 0 0 ; 1 1 0 0 1 0 0 0 0 0 0 1; ...
%     0 0 0 1 0 0 0 0 1 1 1 0 ; 0 1 0 0 0 1 1 0 0 1 0 0; ...
%     1 0 1 0 0 0 0 1 0 0 1 0 ; 0 0 0 1 1 0 0 0 1 0 0 1; ...
%     1 0 0 1 1 0 1 0 0 0 0 0 ; 0 0 0 0 0 1 0 1 0 0 1 1; ...
%     0 1 1 0 0 0 0 0 1 1 0 0];
%Designed =0;


%DESIGNED (comment this section to run previous examples)
Designed =1;
pcm =[1 1 1 1 0 0 1 0 0 0 0 0; 1 0 0 0 1 1 1 1 0 0 0 0;...
     0 1 1 1 0 0 0 1 1 0 0 0; 0 0 1 0 0 1 0 1 0 1 0 1; ...
     0 0 0 1 1 0 0 0 1 1 1 0; 1 0 0 0 0 0 0 0 0 0 1 1];
gm=[1 0 0 0 1 0 1 1 1 0 0 1; 0 0 1 1 1 1 0 0 0 1 1 1;...
    1 1 1 1 0 0 0 1 0 1 0 1; 0 1 0 0 0 1 1 0 1 0 1 1;...
    1 0 1 0 0 1 0 0 1 0 1 0; 0 0 0 1 0 0 1 1 0 0 1 1]; 
%for the last DESIGNED LDPC: e target=0.4 


e = 0.1;

pcm_size = size(pcm);
n = pcm_size(2);
k = n-pcm_size(1);

mistakes_histogram=zeros(1,n+1);
corrected_codewords_histogram =zeros(1,n+1);
perc_corr_histogram=zeros(1,n+1);

%codeword set
codewords_S=codeword_set_f(pcm);

if(Designed ==0)
    gm = gm_from_pcm(pcm);
end

Q = 10000; %number of codewords to send
N = n*Q; %bit stream length

%Transmitter

bit_stream = char(1,N);

for i=1:Q
   
    u_sequence = randi([0 1],k,1)';
    x = mod(u_sequence * gm,2) ;
 
    bit_stream((i-1)*n+1:i*n) = erase(num2str(x),' ');    
    
end

% Channel 
channel_stream = zeros(1,N);
%channel_stream = round((1/(2*(1-e)))*rand(N,1));
for i=1:N
    if(rand<1-e)
        channel_stream(i)=0;
    else
        channel_stream(i)=1;
    end
end

% Receiver
rec_stream = bit_stream;
for i=1:N
    
    if channel_stream(i)==1
        rec_stream(i)='?';
    end
    
end


decoded_stream1 = char(1,N);
decoded_stream2 = char(1,N);
codeword_found_errors=0;
total_codeword_errors=0;
codewords_transmitted_falsely=0;
% Decoding sequence with bit decoder
for i=1:Q
   
   rec_sequence = rec_stream((i-1)*n+1:i*n);
   true_sequence = bit_stream((i-1)*n+1:i*n);
   
   if ~all(true_sequence' == rec_sequence')
       codewords_transmitted_falsely=codewords_transmitted_falsely + 1;
   end
   
   mistakes = numel(rec_sequence(rec_sequence~=true_sequence));
   mistakes_histogram(mistakes+1)=mistakes_histogram(mistakes+1)+1;
    
   %decode
   [decoded_sequence1,cdf]=decode(rec_sequence','01?',[1-e,0,e;0,1-e,e],pcm,[]);
   
   
   if cdf==0
       decoded_sequence1=rec_sequence;
       codeword_found_errors = codeword_found_errors + 1;
   end
   
   
   if ~all(true_sequence' == decoded_sequence1)
       total_codeword_errors=total_codeword_errors+1;
   end
   
   if all(true_sequence' == decoded_sequence1)
         corrected_codewords_histogram(mistakes+1)= corrected_codewords_histogram(mistakes+1)+1;
   end
   
   decoded_stream1((i-1)*n+1:i*n) = decoded_sequence1;
   
   %block map decoding
   %decoded_sequence2 = BlockMapDecoding(rec_sequence', '01?', [1-e,0,e;0,1-e,e],codewords_S);
   
   %decoded_stream2((i-1)*n+1:i*n) = decoded_sequence2;
   
end

for i=1:n+1
        if(mistakes_histogram(i)==0)
            perc_corr_histogram(i)=0;
        else
            perc_corr_histogram(i) = corrected_codewords_histogram(i)/ mistakes_histogram(i);
        end
end


plot_until=8;
plotdata=[100*mistakes_histogram(1:plot_until+1)'/Q, 100*corrected_codewords_histogram(1:plot_until+1)'/Q];
figure
bar(0:plot_until, plotdata, 'grouped')
legend('Percentage of sequences received with q mistakes','Percentage of sequences CORRECTLY retrieved with q mistakes')
title('Percentage histogram of sequences received and retrieved for e=0.1')
xlabel('Number of mistakes q')

figure
stem(0:plot_until,100*perc_corr_histogram(1:plot_until+1))
legend({"Percentage of sequences CORRECTLY retrieved with q mistakes"+newline+"with respect to those that were sent"})
title('"Ability" to decode sequences with q mistakes')
xlabel('Number of mistakes q')
ax=gca;
xlim([-1 10])
ylim([-5 130])
%set(ax, 'YTick', [ 0,20,40,60,80, 100])
ytickformat(ax, 'percentage');


error_stream1=decoded_stream1~=bit_stream;

BER=numel(error_stream1(error_stream1==1))/N; %bit error rate with the belief propagation decoder

