%% Mashup Example
addpath code
clear all
clc
%% Required external dependencies: see README.txt for more information
%addpath /path-to-libsvm-package/matlab % LIBSVM package (for cross_validation.m)
%addpath /path-to-lbfgs-package % L-BFGS package (only if svd_approx = false)

%% Example parameters
org = 'yeast';      % use human or yeast data
onttype = 'level1'; % which type of annotations to use
                    %   options: {bp, mf, cc} for human GO,
                    %            {level1, level2, level3} for yeast MIPS
ontsize = [];       % consider terms in a specific size range (*human GO only*)
                    %   examples: [11 30], [31 100], [101 300]
nperm = 5;          % number of cross-validation trials
svd_approx = true;  % use SVD approximation for Mashup
                    %   recommended: true for human, false for yeast
ndim = 500;         % number of dimensions
                    %   recommended: 800 for human, 500 for yeast

%% Construct network file paths
string_nets = {'neighborhood', 'fusion', 'cooccurence', 'coexpression', ...
               'experimental', 'database'};
network_files = cell(1, length(string_nets));
for i = 1:length(string_nets)
  network_files{i} = sprintf('mashup_data/networks/%s/%s_string_%s_adjacency.txt', ...
                             org, org, string_nets{i});
end

%% Load gene list
gene_file = sprintf('mashup_data/networks/%s/%s_string_genes.txt', org, org);
genes = textread(gene_file, '%s');
ngene = length(genes);

%% Load known annotations
fprintf('[Loading annotations]\n');
if strcmp(org, 'human')
  anno = load_go(onttype, genes, ontsize, true);
elseif strcmp(org, 'yeast')
  anno = load_mips(onttype, genes);
end
fprintf('Number of functional labels: %d\n', size(anno, 1)); 

% Select one categories for regression purpose.
[v,i] = max(sum(anno,2));
annotation_1 = anno(i,:)';

%% Mashup integration
fprintf('[Mashup]\n');
%beta = mashup(network_files, ngene, ndim, svd_approx,anno);
%Concat each network(6400 by 6400 matrix) into a single matrix. 
%I first run the diffusion and then concat each network together.
n_networks = length(network_files);
networks = zeros(ngene * n_networks, ngene);
for i = 1:length(network_files)
    fprintf('Loading %s\n', network_files{i});
    A = load_network(network_files{i}, ngene);%load the similarirty networks.
    fprintf('Running diffusion\n')
    Q = rwr(A, 0.5);%running random walk.
    %R = log(Q + 1/ngene); % smoothing
    start = ngene * (i-1)+1;
    networks(start:(start+ngene-1),:) = Q;%concat each network together.
    %RR_sum = RR_sum + R * R';
end
%%
%because we have 6 netwokrs and 6400 genes
%networks is a (38400-by-6400) matrix, I just concat each diffusion network
%together.
fprintf('Running Integration\n');
clear R Q A
%Then, we run the pca to obtain eigenvectors and eigenvalue. Internally,
%this function would first de-mean the networks matrix, 
%compute the covariance matrix and then perform SVD algorithm on it.
%eigenvalue is actually stored in a decreasing order in latent, 
%so dose the coeff matrix(each column represents a eigenvector).
%Here, coeff is the 6400x6400 matrix, score is the networks matrix mapped
%into eigen-space(From your description, we do not need this), and latent
%is the 6400 dimensional eigenvalues vector.
tic
[coeff,score,latent] = pca(networks,'Algorithm','eig');
save('result/mashup/eigenvalue.mat','latent')
fprintf('PCA finished, computing beta vector')
toc
%annotation_1 is 6400 dimensianl vector(binary value for
%regression).coeff is the 6400x6400 eigenvector matrix.
%Here, we want to run a linear regression, that's coeff X beta = annotation_1.
%I simply do a left divide to obtain beta vector (6400 dimension).
tic
beta = coeff \ annotation_1;
toc
%Finally, we have both the beta and eigenvalue(latent variable). 
%We need to compare the magtitude of beta and eigenvalue, so I do some
%normalization to make the comparison more intuitive.
%It must be mentioned that beta contains a lot of negative value while
%eigenvalue is all positive. So I simply take the absolute value of the
%beta.

save('result/mashup/linear_regression_input.mat','coeff','annotation_1')
save('result/mashup/linear_regression_output.mat','beta')
beta_ = abs(beta)/sum(abs(beta))/7;%7 is just a number choosen to scale for plot.

latent_ = latent/sum(latent);
%% Plot
plot(1:ngene,beta_,1:ngene,latent_)
legend('beta','eigenvalue')
saveas(gcf,'result/mashup/comparison.png')
