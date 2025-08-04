function out=num_codeword(in)

    len=length(in);
    
    intermidiate=char(1,2*len-1);
    
    for i=1:len
        intermidiate(2*i-1)=in(i);
    end

    for i=1:len-1
        intermidiate(2*i)=' ';
    end

    out=str2num(intermidiate)';

end