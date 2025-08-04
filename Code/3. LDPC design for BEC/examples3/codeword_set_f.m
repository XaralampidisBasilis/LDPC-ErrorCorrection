function out_set=codeword_set_f(pcm)

    pcm_size = size(pcm);
    bits_number = pcm_size(2);
    
    index=1;
    
    out_set=zeros(2^bits_number-1,bits_number);
    
    for i = 0 : 2 ^ bits_number - 1

        codeword = num_codeword(dec2bin(i, bits_number));
        result = ~any(mod(pcm * codeword', 2));
        
        if result
            
            out_set(index,:)=codeword;
            index =index+1;
          
        end
        
    end

    out_set=out_set(1:index-1,:);


end





