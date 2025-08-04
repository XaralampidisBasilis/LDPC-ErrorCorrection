function transmited_codeword = BlockMapDecoding(received_codeword, receiver_symbols, Channel_Matrix, ...
codeword_set)
%DECODE3 Maximum Likelihood Criterion for Block Map Decoding 
%                               Input arguments
%--------------------------------------------------------------------------------
%{                        
   receiver_symbols = a char array contaning the symbols that the receiver
                     has in proper order as the Channel_Matrix
                     ex: y in ['0' '1' '?'] = '01?'

   Channel_Matrix = a matrix containing the dependet probabilities of the
                    transmition channel, Channel_Matrix(i,j) = PY|X(yj|xi)
                    ex:   xy  0    1    ?    
                          0   1-å  0    å
                          1   0    1-å  å

   Bipartite_Matrix = contains information about the bipartime graph 
                      which relates the factorization functions with their
                      respective arguments. This logic matrix has ones
                      when the functions located in the rows relate to
                      their arguments located in the columns
                      
                      ex:       x1   x2
                             f1  1    1
                             f2  0    1 

   received_codeword = the codeword that the receiver detects with errors
                       Y = (y1 y2 y3 ... yn)
                       ex: Y = '01?101?010'

%}
%                               Output arguments
%--------------------------------------------------------------------------------
%{
   transmited_codeword = the information codeword that the trasmier sends 
                          through the channel X = (x1 x2 x3 ... xn)
                         ex: X = '0101011010'

%}
    

    symbols_index = sum((received_codeword == receiver_symbols') .* (1:length(receiver_symbols))', 1);
    P = Channel_Matrix(:, symbols_index);
    
    bits_number = length(received_codeword);
    codewords_number = length(codeword_set);
    products_array = zeros(codewords_number, 1);

    for i = 1 : codewords_number
        
            codeword=codeword_set(i,:);
            
            product = P(codeword(1) + 1, 1);
            for b = 2 : bits_number
                product = product * P(codeword(b) + 1, b);
            end
            products_array(i ) = product;
        
    end

    [~, index] = max(products_array);
    transmited_codeword = sprintf('%d', codeword_set(index,:));
end

