%Code that shows the difference in calculation complexity when
%using factor graphs compared to not

N=30; 
% warning N increases complexity exponentially N=7 ~ 5 seconds
% -fixed: the brute force algorithm caluclates only for N<12


bm=[1 0 1 1 0 0; 1 0 0 0 0 0; 1 1 0 0 1 0; 0 0 1 0 0 1]; 
variable_name='x4';
f1=@(x1,x3,x4) (x1+x3)/((x4-x1).^2+1);
f2=@(x1) x1.^2;
f3=@(x1,x2,x5) x5.^x2+x1;
f4=@(x3,x6) x3+x6.^2+1;
f={f1,f2,f3,f4};
F={@(x1,x2,x3,x4,x5,x6) f1(x1,x3,x4)*f2(x1)*f3(x1,x2,x5)*f4(x3,x6)};
B=[1 1 1 1 1 1];
hold off

times = zeros(2,N-1);


for i =2:N
    
    variable_vector=(1:i);
    
    tic
    factor_graph_marginal = marginal(variable_vector,variable_name,bm,f, @summation, @multiplication);
    times(1,i-1)=toc;
    
    if(i<12)
    tic
    classic_marginal = marginal(variable_vector,variable_name,B,F, @summation, @multiplication);
    times(2,i-1)=toc;
    end
    
end

%plotting
t1=2:N;
t2=2:min(N,11);
figure
hold on
plot(t1,times(1,:))
plot(t2,times(2,1:min(N,10)))
set(gca, 'YScale', 'log')
legend('Factor graph algorithm','Brute-force algorithm')
title('Algorithm time comparison')
xlabel('Set cardinality |X|')
ylabel('Calculation time (s)')