function [AAA] = FOS_LTI(A,alpha,N,n,p)
    
    N = p;  % The interval between N and p is not used for AAA, reduce computation time

    %Set up FOS LTI A matrix
    % Setting up A_0, A_1, ...
    AA = cell(N+1,1);
    D = zeros(n,n);
    j = 0;
    for i=1:n
        D(i,i) = 1;
        for l=0:j
            D(i,i) = D(i,i)*(j-l-alpha(i))/(l+1);
        end
    end
    AA{1} = A - D;

    for j=1:(N+1)
        D = zeros(n,n);
        for i=1:n
            D(i,i) = 1;
            for l=0:j
                D(i,i) = D(i,i)*(j-l-alpha(i))/(l+1);
            end
        end
        AA{j+1} = -D;
    end

    % Setting up parameters in augmented LTI system (i.e. finite-history FOS)
    AAA = [];
    for k=1:p
        AAA = [AAA AA{k}];
    end
    for j=2:p
        AAAj = [];
        for k=1:p
            if j==k+1
                AAAj = [AAAj eye(n)];
            else
                AAAj = [AAAj zeros(n,n)];
            end
        end
        AAA = [AAA; AAAj];
    end

end

