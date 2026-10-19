function H = haltonseq(N, s)
    % Generates N Halton points of dimension s
    p = haltonset(s);
    H = net(p, N);
end
