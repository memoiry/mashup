
type PatientsNetworks

    index

end

function load_net_file_name(string_nets, dir=nothing)::Array{String,1}
    if dir != nothing
        for i = 1:length(string_nets)
            string_nets[i] = joinpath(dir,string_nets[i])
        end
    end
end


function load_patients_index(filename::String,
                             )::Dict{String,Int32}



end


function load_disease(filename::Stirng,
                      )::Matrix


end

dump(eigenvalue_list::Matrix, β::Vector) =  

```
Mashup algorithm
```

function pca(A::Matrix, num::Int64)
    n = size(A,1)
    A_demean = broadcast(/,A,mean(A,1))
    Σ = 1 / n * A_demean' * A_demean
    eigenvalue, eigenvector =  eig(Σ)
    return eigenvalue, eigenvector
end

function load_net(filename::String,
                  database::GMANIA)
    network = readdlm(filename);
    n_relation = size(network, 1)
    if database.use_index
        for i in 1:n_relation
            for j in 1:2
                network[i,j] = database.patients_index[network[i,j]] 
            end
        end
    end
    if isa(network, Array{any,2}) 
        network = convert(Array{Float64,2}, network)
    end
    @assert isa(network, Array{Float64, 2]})
    A = full(sparse(M(:,1), M(:,2), M(:,3), database.n_patiens, database.n_patiens));
    if !isequal(A, A') # symmetrize
        A = A + A'
    end
    A = A + diagm(sum(A, 2) .== 0)
end

function mashup(net_file_name::Array{String,1},
                n_patiens::Int32,
                            )
    n_net = length(net_file_name)
    net = zeros(n_patiens * n_net, n_patiens);
    eigen_value_list = zeros(n_patiens,length(net_files));
    for i = 1:length(net_files)
        fprintf('Loading %s\n', net_files{i});
        A = load_net(net_files{i}, n_patiens, index_table);#load the similarirty net.
        fprintf('Running diffusion\n')
        Q = rwr(A, 0.5);#running random walk.
        #R = log(Q + 1/n_patiens); % smoothing
        start = n_patiens * (i-1)+1;
        net(start:(start+n_patiens-1),:) = Q;%concat each net together.
        [coeff,score,latent] = pca(A,'Algorithm','eig');
        if length(latent) == n_patiens - 1
            latent = [latent;0];
        end
    eigen_value_list(i,:) = latent;
    #RR_sum = RR_sum + R * R';
end

    return eigenvalue_list, β
end

```
Random walk with restart
```

function rwr(A::Matrix, restart_prob)
    n = size(A, 1);
    A = A - diagm(diag(A));
    A = A + diagm(sum(A,1) .== 0);
    P = broadcast(/, A, sum(A,1))
    Q = (eye(n) - (1 - restart_prob) * P) \ (restart_prob * eye(n));
end


function mashup_test()
    netwokrs_dir = 
    net_file_name = load_net_file_name()
    patients_index, n_patients = load_patients_index()
    disease = load_disease()
    eigenvalue_list, β = mashup(net_file_name, patients_index, disease)

end

