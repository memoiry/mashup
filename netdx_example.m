%% netDX example
clear all
clc
addpath code

%% Construct network file paths
string_nets = {'age', 'grade', 'stage'};
network_files = cell(1, length(string_nets));
for i = 1:length(string_nets)
    network_files{i} = sprintf('netDX_data/KIRC/data/networks/%s.txt', string_nets{1});
end
%% Construct index table to convert the patient name to index key.
%Here we have 150 patients.
fid = fopen('netDX_data/KIRC/KIRC_clinNets_pheno_matrix.txt');
dataset = textscan(fid,'%d%s')
genes_name = dataset{2};
genes_index = dataset{1};
ngene = length(dataset{1});
index_table = containers.Map;
for i = 1:ngene
    index_table(genes_name{i}) = genes_index(i);
end
%% Load known annotations

anno = dlmread('netDX_data/KIRC/data/annotations/networks.csv');
annotation = anno(:,2);

%% Mashup integration
fprintf('[Mashup]\n');
%beta = mashup(network_files, ngene, ndim, svd_approx, anno');
n_networks = length(network_files);
networks = zeros(ngene * n_networks, ngene);
for i = 1:length(network_files)
    fprintf('Loading %s\n', network_files{i});
    A = load_network_netdx(network_files{i}, ngene, index_table);%load the similarirty networks.
    fprintf('Running diffusion\n')
    Q = rwr(A, 0.5);%running random walk.
    %R = log(Q + 1/ngene); % smoothing
    start = ngene * (i-1)+1;
    networks(start:(start+ngene-1),:) = Q;%concat each network together.
    %RR_sum = RR_sum + R * R';
end
fprintf('Running Integration\n');
clear R Q A
[coeff,score,latent] = pca(networks,'Algorithm','eig');
save('result/netdx/eigenvalue.mat','latent')
fprintf('PCA finished, computing beta vector')
beta = coeff \ annotation;
save('result/netdx/linear_regression_input.mat','coeff','annotation')
save('result/netdx/linear_regression_output.mat','beta')
beta_ = abs(beta)/sum(abs(beta))/2;
latent_ = latent/sum(latent);
%% Plot
plot(1:ngene,beta_,1:ngene,latent_)
legend('beta','eigenvalue')
saveas(gcf,'result/netdx/comparison.png')
