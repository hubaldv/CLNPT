function [AAAA,BBBB]=predmodgen(AAA,BBB,P)

%Prediction matrices generation
%This function computes the prediction matrices to be used in the
%optimization problem

% Setting up auxiliary matrices to compute QP solution
    AAAA = [];
    for k=1:P
        AAAAj = [];
        for j=1:P
            if k<j
                AAAAjk = eye(size(AAA,1));
                for l=k:(j-1)
                    AAAAjk = AAA*AAAAjk;
                end
            elseif k==j
                AAAAjk = eye(size(AAA,1));
            else
                AAAAjk = zeros(size(AAA,1),size(AAA,1));
            end
            AAAAj = [AAAAj; AAAAjk];
        end
        AAAA = [AAAA AAAAj];
    end

    BBBB = [];
    for j=1:P
        BBBBj = [];
        for k=1:P
            if j==k
                BBBBjk = BBB;
            else
                BBBBjk = zeros(size(BBB,1),size(BBB,2));
            end
            BBBBj = [BBBBj; BBBBjk];
        end
        BBBB = [BBBB BBBBj];
    end

end


