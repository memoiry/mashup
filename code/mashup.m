%% network_files (cell array): paths to adjacency list files
%% ngene (int): number of genes in input networks
%% ndim (int): number of output dimensions
%% svd_approx (bool): whether to use SVD approximation for large-scale networks
%%
function beta = mashup(network_files, ngene, ndim, svd_approx, anno)
if svd_approx
    n_networks = length(network_files);
    networks = zeros(ngene * n_networks, ngene);
    for i = 1:length(network_files)
        fprintf('Loading %s\n', network_files{i});
        A = load_network(network_files{i}, ngene);
        fprintf('Running diffusion\n')
        Q = rwr(A, 0.5);
        
        %R = log(Q + 1/ngene); % smoothing
        start = ngene * (i-1)+1;
        networks(start:(start+ngene-1),:) = Q;
        %RR_sum = RR_sum + R * R';
    end
    
    fprintf('Running Integration\n');
    clear R Q A
    tic
    [coeff,score,latent] = pca(networks,'Algorithm','eig');
    toc
    fprintf('PCA finished, computing beta vector')
    [v,i] = max(sum(anno,2));
    annotation_1 = anno(i,:)';
    tic
    %size(annotation_1)
    %size(networks)
    beta = coeff \ annotation_1;
    beta_ = abs(beta)/sum(abs(beta))/7;
    latent_ = latent/sum(latent);
    plot(1:ngene,beta_,1:ngene,latent_)
    
    toc
    
    %fprintf('All networks loaded. Learning vectors via SVD...\n');
    %[V, d] = eigs(cov_, ndim);
    %x = diag(sqrt(sqrt(diag(d)))) * V';
else
    Q_concat = [];
    for i = 1:length(network_files)
        fprintf('Loading %s\n', network_files{i});
        A = load_network(network_files{i}, ngene);
        fprintf('Running diffusion\n');
        Q = rwr(A, 0.5);
        
        Q_concat = [Q_concat; Q];
    end
    clear Q A
    Q_concat = Q_concat / length(network_files);
    
    fprintf('All networks loaded. Learning vectors via iterative optimization...\n');
    x = vector_embedding(Q_concat, ndim, 1000);
end

fprintf('Mashup features obtained.\n');
    function A = load_network(filename, ngene)
        M = dlmread(filename);
        A = full(sparse(M(:,1), M(:,2), M(:,3), ngene, ngene));
        if ~isequal(A, A') % symmetrize
            A = A + A';
        end
        A = A + diag(sum(A, 2) == 0);
    end
end
