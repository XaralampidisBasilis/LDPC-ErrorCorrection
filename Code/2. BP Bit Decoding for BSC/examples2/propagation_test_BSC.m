


pcm = [ 1 1 1 1 0 0; 0 1 1 0 1 0; 1 0 1 0 0 1];
% pcm=[0 0 1 0 0 1 1 1 0 0 0 0 ; 1 1 0 0 1 0 0 0 0 0 0 1; ...
%     0 0 0 1 0 0 0 0 1 1 1 0 ; 0 1 0 0 0 1 1 0 0 1 0 0; ...
%     1 0 1 0 0 0 0 1 0 0 1 0 ; 0 0 0 1 1 0 0 0 1 0 0 1; ...
%     1 0 0 1 1 0 1 0 0 0 0 0 ; 0 0 0 0 0 1 0 1 0 0 1 1; ...
%     0 1 1 0 0 0 0 0 1 1 0 0];

e = 0.02;

pcm_size = size(pcm);
n = pcm_size(2);
k = n-pcm_size(1);

gm = gm_from_pcm(pcm);

Q = 5000; %number of codewords to send
N = n*Q; %bit stream length

bit_stream = zeros(N,1);

for i=1:Q
   
    u_sequence = randi([0 1],k,1)';
    x = mod(u_sequence * gm,2) ;
 
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

decoded_stream =zeros(N,1);
codeword_found_errors=0;
total_codeword_errors=0;
% Decoding sequence with bit decoder
for i=1:Q
   
   rec_sequence = rec_stream((i-1)*n+1:i*n);
   true_sequence = bit_stream((i-1)*n+1:i*n);
   
   pm = pm_from_rec(rec_sequence,e); 
   
   [decoded_sequence,err]=bitDecode(pcm,pm,15);
   if err==1
       decoded_sequence=rec_sequence;
       codeword_found_errors = codeword_found_errors + 1;
   end
   
   if ~all(true_sequence == decoded_sequence)
       total_codeword_errors=total_codeword_errors+1;
   end
   
   decoded_stream((i-1)*n+1:i*n) = decoded_sequence;
   
end

error_stream=mod(decoded_stream-bit_stream,2);
error_wod_stream=mod(rec_stream-bit_stream,2);

BER=numel(error_stream(error_stream==1))/N; %bit error rate with the belief propagation decoder
BER_WOD=numel(error_wod_stream(error_wod_stream==1))/N; %bit error rate WITHOUT the belief propagation decoder(received=decoded)


