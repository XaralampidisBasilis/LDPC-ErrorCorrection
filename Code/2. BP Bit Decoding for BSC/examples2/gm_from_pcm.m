function output_gm=gm_from_pcm(input_pcm)
    % This function outputs a generator matrix corresponding to a given
    % parity check matrix. The algorithms used are taken from Raymond Hill, 
    % A First Course In Coding Theory, OUP, 1986 (pages 50-51) and some
    % observations("Construction of Generator Matrix for a given Parity
    % Check Matrix") that were sent to the students by prof Leonidas
    % Georgiadis.This part of the code may contain mistakes but was a bonus
    % part of the assignment.
    
    pcm_size = size(input_pcm);
 
    pcm_rows = pcm_size(1);
    pcm_cols = pcm_size(2);
    n = pcm_cols;
    k = pcm_cols-pcm_rows;

    norm_pcm = input_pcm;
    
    L = eye(pcm_rows);
    R = eye(pcm_cols);
    Li = eye(pcm_rows);
    Ri = eye(pcm_cols);

    %obtaining the normalized pcm while saving the left and right linear
    %operations we are using
    for j=n:-1:k+1
        
        flag_gjj_is0 = (norm_pcm(j-k,j)==0);
        rj=j-k;
        
        %step 1
        if (flag_gjj_is0)
            i=j-k-1;
            flag_row_found=0;
            if i>0
                while(norm_pcm(i,rj)==0 & i>1)
                    i=i-1;
                end
                if norm_pcm(i,rj)==1
                   flag_row_found=1;
                end
            end
            if flag_row_found
                [norm_pcm(i,:),norm_pcm(rj,:)]=deal(norm_pcm(rj,:),norm_pcm(i,:));
                [L(i,:),L(rj,:)]=deal(L(rj,:),L(i,:));
            else
                h=j-1;
                cj=j;
                flag_col_found=0;
                while(norm_pcm(rj,h)==0 & i>1)
                    i=i-1;
                end
                if norm_pcm(rj,h)==1
                   flag_col_found=1;
                end
                if flag_col_found
                    [norm_pcm(:,cj),norm_pcm(:,h)]=deal(norm_pcm(:,h),norm_pcm(:,cj));
                    [R(h,:),R(cj,:)]=deal(R(cj,:),R(h,:));
                else
                    disp('Parity check matrix not designed correctly.');
                end
            end
        end
        %step 2 ommitted because we are in Field wih base 2
        
        %step 3
        Li = eye(pcm_rows);
        for i=1:pcm_rows
            if(i~=rj)
                if(norm_pcm(i,j)==1)
                    Li(i,rj)=1;
                end
                norm_pcm(i,:)=mod(norm_pcm(i,:)-norm_pcm(i,j)*norm_pcm(rj,:),2);
            end
        end
        L=mod(Li*L,2);
    end

    %constructing the generator matrix
    
    B=norm_pcm(:,1:k);
    
    norm_gm=[eye(k) B'];

    output_gm = mod(norm_gm*R',2);

end