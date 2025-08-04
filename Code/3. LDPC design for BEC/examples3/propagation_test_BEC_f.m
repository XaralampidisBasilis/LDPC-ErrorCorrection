function [BERS, mistakes_histogram,corrected_codewords_histogram] = propagation_test_BEC_f(pcm,gm,e,Q)

    BERS=zeros(1,3);
    %BERS 1 decode full errors
    %BERS 2 decode "guessed wrong bit" errors
    %BERS 3 Block Map Decoding errors
  
    
    
    pcm_size = size(pcm);
    n = pcm_size(2);
    k = n-pcm_size(1);
    
    mistakes_histogram=zeros(1,n+1);
    corrected_codewords_histogram =zeros(1,n+1);
    perc_corr_histogram=zeros(1,n+1);
    
    
    
    %codeword set
    codewords_S=codeword_set_f(pcm);
    
    
    
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
    % Decoding sequence with bit decoder
    for i=1:Q
   
        rec_sequence = rec_stream((i-1)*n+1:i*n);
        true_sequence = bit_stream((i-1)*n+1:i*n);
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
        
        %Block map decoding
         
        decoded_sequence2 = BlockMapDecoding(rec_sequence', '01?', [1-e,0,e;0,1-e,e],codewords_S);
        decoded_stream2((i-1)*n+1:i*n) = decoded_sequence2;
        
    end

    error_stream1=decoded_stream1~=bit_stream;
    error_stream2=decoded_stream2~=bit_stream;
    
    BERS(1)=numel(error_stream1(error_stream1==1))/N;
    BERS(2)=numel(error_stream1(error_stream1==1 & decoded_stream1~='?'))/N;
    BERS(3)=numel(error_stream2(error_stream2==1))/N;

    for i=1:n+1
        if(mistakes_histogram(i)==0)
            perc_corr_histogram(i)=0;
        else
            perc_corr_histogram(i) = corrected_codewords_histogram(i)/ mistakes_histogram(i);
        end
    end
end