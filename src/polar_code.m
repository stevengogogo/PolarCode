function polar_code
l = 1;
e = 0.5;
n = 2^l;
T_n = Transform(l)
I_SUM = MutualInformation(T_n, e)
end

function M = Transform(l)
    M = Permut(l)*Transfer(l);
end

function p = Permut(l)
    p = zeros(2^l);
    for i = 1:2^l
        p(index_naming(i,l),i) = 1;
    end
end

function T = Transfer(l)
    T = [1,1;1,0];
    for i = 2:l
        T = kron(T,T);   
    end
end

function index = index_naming(i,l)
    index = bi2de(fliplr(de2bi((i-1),l)))+1;
end

function I = MutualInfo_vec(w_vec,e)
    n = nnz(w_vec);
    l = length(w_vec);
    I = 0;
    for i = 0:n
        I = I + nchoosek(l,i)*e^i*(1-e)^(l-i); 
    end
end

function I_sum = MutualInformation(T_n, e)
    I_sum = 0;
    for i = 1:length(T_n)
        I_sum = I_sum + MutualInfo_vec(T_n(i,:),e);
    end
end