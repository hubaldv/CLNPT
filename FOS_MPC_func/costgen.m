function [Htotal,ftotal, QQ, cc]=costgen(lambda,lambdaR,AAA,BBB,AAAA,BBBB,xxx0,n,P)

    % Q, R, c matrices/vectors
    Q = lambda*zeros(size(AAA,1));
    Q(1:n,1:n) = lambda*eye(n);
    R = lambdaR*eye(size(BBB,2));
    c = zeros(size(AAA,1),1);

    QQ = [];
    for j=1:P
        Qj = [];
        for k=1:P
            if j==k
                Qjk = Q;
                % uncomment to induce seizure
                %if j~=1
                %   Qjk = 0*Qjk;
                %end
            else
                Qjk = zeros(size(Q,1),size(Q,1));
            end
            Qj = [Qj; Qjk];
        end
        QQ = [QQ Qj];
    end


    RR = [];
    for j=1:P
        Rj = [];
        for k=1:P
            if j==k
                Rjk = R;
            else
                Rjk = zeros(size(BBB,2),size(BBB,2));
            end
            Rj = [Rj; Rjk];
        end
        RR = [RR Rj];
    end

    cc = [];
    for k=1:P
        cc = [cc; c];
    end

    Htotal = 2*(BBBB'*AAAA'*QQ*AAAA*BBBB + RR);
    ftotal = 2*BBBB'*AAAA'*QQ*AAAA(:,1:size(AAA,1))*AAA*xxx0 + BBBB'*AAAA'*cc;
end