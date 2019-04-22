function B = normaleq(X,Y)

    B = pinv(X' * X) * (X' * Y);
    
end