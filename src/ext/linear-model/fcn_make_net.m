function A = fcn_make_net(omega,ci,M,N)
count = 0;
A = zeros(N);
I = dummyvar(ci);
while count < M
    r = randi(K);
    s = randi(K);
    d = I(:,r)*I(:,s)';
    d = d | d';
    if rand < omega(r,s)
        [u,v] = find(triu(d & ~A,1));
        if ~isempty(u)
            t = randi(length(u));
            A(u(t),v(t)) = 1;
            count = count + 1;
        end
    end
end
A = A + A';