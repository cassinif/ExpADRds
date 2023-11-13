function x = thom_ad(ldiag,ddiag,udiag,w,b)
    % ldiag has initial leading zeros
    % udiag has final leading zeros

    n = length(ldiag);
    alphai = NaN(n,1);
    betai = NaN(n,1);
    x = NaN(n,1);

    for i=1:w
        alphai(i) = udiag(i)/ddiag(i);
        betai(i) = b(i)/ddiag(i);
    end

    for i=(w+1):(n-w)
        alphai(i) = udiag(i)/(ddiag(i)-alphai(i-w)*ldiag(i));
    end

    for i=(w+1):n
        betai(i) = (b(i)-betai(i-w)*ldiag(i))/(ddiag(i)-alphai(i-w)*ldiag(i));
    end

    for i=1:w
      x(n-i+1) = betai(n-i+1);
    end

    for i=(n-w):-1:1
        x(i) = betai(i) - alphai(i)*x(i + w);
    end
end
