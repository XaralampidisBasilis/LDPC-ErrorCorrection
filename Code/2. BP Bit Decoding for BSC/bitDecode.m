function [decoded_sequence ,error]=bitDecode(paritycheck_matrix,probabilities_matrix,N_max)
    %This function models an algorithm of bit-decoding for error correction. 
    %It takes 3 inputs:
    %   -paritycheck_matrix: Is the parity check matrix of the given code.
    %   We could also imagine that it is the bipartite matrix of the
    %   function-node graph(the bipartite matrix has 1 when a variable node is
    %   connected to a fucntion node).
    %   -probabilities_matrix: Is the matrix that corresponds to
    %   [p(y=y_known|x=0);p(y=y_known|x=1)] which implements channel-state
    %   information(y_known is the bit that the receiver believes was sent in the first place). 
    %   In reality we only use the log of the ratio of the two probabilities(log likelihood ratio).
    %   -N_max: Is the number of iterations that the algorthm implements, before it declares an error.
    %   
    %It produces 2 ouputs:
    %   -decoded_sequence: The most probable sended sequence according to
    %   the algorithm
    %   -error: a bit that is 0 when the algorithm has found a codeword and
    %   1 if not, hence an error occured
    %
    %A quick view of the algorithm is:
    %   1.Find and entitle the edges in the bipartite graph(to lower 
    %   the computing cost of the algorithm).
    %   2.Make two maps for variable and function messages:check_msg_needed
    %   and var_msg_needed.The maps determine which edges should be computed
    %   to determine the output in each edge(coming from variables or functions 
    %   respectively) in the main iteration.  
    %   3.Implement the N_max message passing iteration between function and
    %   variable nodes.
    %
    %  A possible call is:
    %  pcm=[1,1,1,1,0,0;0,0,1,1,0,1;1,0,0,1,1,0];
    %  pm=[0.01 0.99; 0.02 0.98; 0.03 0.97; 0.6 0.4; 0.95 0.05; 0.94 0.06];
    %  [dec_seq , er]= bitDecode(pcm,pm,10);


    %1.We encapsulate the graph information(which variable node is connected to
    %which check node) into a matrix called edges.The first part corresponds to
    %the check node number and the second one to the variable node number
    pcm_size = size(paritycheck_matrix);
    function_No = pcm_size(1);
    variable_No = pcm_size(2);

    edges_No = sum(paritycheck_matrix(:)==1);
    edges = zeros(edges_No,2);

    u=1;%u is an indice that helps us determine the edges matrix

    for i=1:function_No
        for j=1:variable_No
            if(paritycheck_matrix(i,j)==1)
                edges(u,1)=i;
                edges(u,2)=j;
                u=u+1;
            end
        end
    end

    funcs=edges(:,1);
    vars=edges(:,2);

    %2.Creating a cell structure that stores the values of the edges, messages from which must be 
    %computed in order to determine the output, for each edge.

    var_msg_needed=cell(edges_No,1);
    for i=1:edges_No
        f_node=funcs(i);
        v_node=vars(i);

        indice = find(funcs-f_node==0);
        var_msg_needed{i} = indice(vars(indice)~=v_node);    
    end

    check_msg_needed=cell(edges_No,1);
    for i=1:edges_No
        f_node=funcs(i);
        v_node=vars(i);

        indice = find(vars-v_node==0);
        check_msg_needed{i} = indice(funcs(indice)~=f_node);    
    end

    cmn_for_variable=cell(variable_No,1);%check message needed for variable
    for i=1:variable_No
        v_node=i;

        cmn_for_variable{i} = find(vars-v_node==0);
    end


    %Decleration of matrices needed for the message passing algorithm
    
    decoded=zeros(variable_No,1);
    output=zeros(variable_No,1);

    channel_info = zeros(variable_No,1);% initial messages that stand for channel information
    var_messages = zeros(edges_No,1);% messages that are transferred from variable to check nodes
    fun_messages = zeros(edges_No,1);% messages that are transferred from check to variable nodes

    for i=1:variable_No
        channel_info(i) = log(probabilities_matrix(i,1)/probabilities_matrix(i,2));
    end

    var_messages = channel_info(vars);% initalization of variable messages


    %3.MESSAGE PASSING ITERATION

    for n=1:N_max
        
        %check nodes
        for i=1:edges_No
            fun_messages(i) = 2*atanh(prod(tanh(var_messages(var_msg_needed{i})/2)));
        end


        %variable nodes
        
        %a.decode sequence according to messages from ALL nodes
        for i=1:variable_No
            output(i) = sum(fun_messages(cmn_for_variable{i})) + channel_info(i);
            if(output(i)>=0)
                decoded(i)=0;
            else
                decoded(i)=1;
            end
        end
        
        %b.if decoded sequence is a codeword STOP
        if(~any(mod(paritycheck_matrix*decoded,2)))
            decoded_sequence=decoded;
            error=0;
            return
        end

        %c.determine output messages to check nodes in all edges
        for i=1:edges_No
                var_messages(i) = sum(fun_messages(check_msg_needed{i}))+channel_info(vars(i));
        end


    end

        %if iteration is completed without founding a codeword declare
        %error
        decoded_sequence=decoded;
        error=1;

end