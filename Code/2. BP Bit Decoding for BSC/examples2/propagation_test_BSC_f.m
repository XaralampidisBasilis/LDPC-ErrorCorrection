function [BERS,times]= propagation_test_BSC_f(pcm,e,Q)
    
times=zeros(1,3);
BERS=zeros(1,3);

pcm_size = size(pcm);
n = pcm_size(2);
k = n-pcm_size(1);


gm = gm_from_pcm(pcm);

N = n*Q; %bit stream length


bit_stream = zeros(N,1);

for i=1:Q
   
    u_sequence = randi([0 1],k,1)';
    x = mod(u_sequence * gm,2) ;
    
    while(any(mod(pcm*x',2)))
        u_sequence = randi([0 1],k,1)';
        x = mod(u_sequence * gm,2) ;
    end
 
    bit_stream((i-1)*n+1:i*n) = x;    
    
end

% Channel 
channel_stream = zeros(N,1);
%channel_stream = round((1/(2*(1-e)))*rand(N,1));
for i=1:N
    if(rand<1-e)
        channel_stream(i)=0;
    else
        channel_stream(i)=1;
    end
end


% Receiver
rec_stream = mod(bit_stream+channel_stream,2);

decoded_stream1 =zeros(N,1);
codeword_found_errors1=0;

decoded_stream2 =zeros(N,1);
codeword_found_errors2=0;

decoded_stream3 =zeros(N,1);





% Decoding sequence with decoder
for i=1:Q
   
   rec_sequence = rec_stream((i-1)*n+1:i*n);
   pm = pm_from_rec(rec_sequence,e); 
   
   %bitDecode
   tic
   [decoded_sequence1,err1]=bitDecode(pcm,pm,15);
   times(1)=times(1)+toc;
   if err1==1
       decoded_sequence1=rec_sequence;
       codeword_found_errors1 = codeword_found_errors1 + 1;
   end
   
   decoded_stream1((i-1)*n+1:i*n) = decoded_sequence1;
   
   %decode
   tic
   [decoded_sequence2,cdf2]=decode(erase(num2str(rec_sequence)',' '),'01',[1-e,e;e,1-e],pcm,[]);
   times(2)=times(2)+toc;
   decoded_sequence2=num_codeword(decoded_sequence2);
   if cdf2==0
       decoded_sequence2=rec_sequence;
       codeword_found_errors2 = codeword_found_errors2 + 1;
   end
   
   decoded_stream2((i-1)*n+1:i*n) = decoded_sequence2;
   
   
   %Block Map Decoding
   tic
   decoded_sequence3=BlockMapDecoding(erase(num2str(rec_sequence)',' '),'01',[1-e,e;e,1-e],pcm);
   times(3)=times(3)+toc;
   decoded_sequence3=num_codeword(decoded_sequence3);
   
   decoded_stream3((i-1)*n+1:i*n) = decoded_sequence3;
   
end

    times=times/Q;
    
    error_stream1=mod(decoded_stream1-bit_stream,2);
    error_stream2=mod(decoded_stream2-bit_stream,2);
    error_stream3=mod(decoded_stream3-bit_stream,2);
    
    
    BERS(1)=numel(error_stream1(error_stream1==1))/N; %bit error rate
    BERS(2)=numel(error_stream2(error_stream2==1))/N;
    BERS(3)=numel(error_stream3(error_stream3==1))/N;
    





end