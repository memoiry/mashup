    function A = load_network_netdx(filename, ngene, index_table)
        fid = fopen(filename);
        network = textscan(fid, '%s%s%f');
        n_relation = length(network{1});
        M = zeros(n_relation,3);
        for i = 1:n_relation
            M(i,1) = index_table(network{1}{i});
            M(i,2) = index_table(network{2}{i});
            M(i,3) = network{3}(i);
        end
        A = full(sparse(M(:,1), M(:,2), M(:,3), ngene, ngene));
        if ~isequal(A, A') % symmetrize
            A = A + A';
        end
        A = A + diag(sum(A, 2) == 0);
    end