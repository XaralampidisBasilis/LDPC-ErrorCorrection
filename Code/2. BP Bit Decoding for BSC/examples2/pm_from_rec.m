function probabilities_matrix=pm_from_rec(received_codeword,e)
    
    %This function gives the probability matrix of p(x=0 or 1 |y=received) 0 or 1 in
    %each column, given that our message goes through a BSC and the
    %receiver knows that.

    len=length(received_codeword);
    probabilities_matrix=zeros(len,2);
    
    for i=1:len
        
        if received_codeword(i)==0
            probabilities_matrix(i,1)=1-e;
            probabilities_matrix(i,2)=e;
        else
            probabilities_matrix(i,1)=e;
            probabilities_matrix(i,2)=1-e;
        end
        
    end

end