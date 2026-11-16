function profile=critical_mode_degree_profile(layers,L_inter,ratio,lambda_star)
% Cumulative critical-mode mass, ranked by mean intra-layer node degree.
% Repeated eigenvalues use the normalized eigenspace projector, so the
% result is independent of the particular orthonormal eigenbasis.
N=size(layers{1},1); K=numel(layers);
A=sparse(blkdiag(layers{:})+ratio*kron(sparse(L_inter),speye(N)));
tolerance=1e-8*max(1,abs(lambda_star));
% The c-panel has at most 2500 layer nodes with the default parameters.
% A complete decomposition also captures every vector of a repeated mode.
[V,D]=eig(full(A));
chosen=abs(diag(D)-lambda_star)<=tolerance;
V=V(:,chosen);
assert(~isempty(V),'Cannot recover the cached critical eigenvalue.');
[V,~]=qr(V,0);
residual=norm(A*V-lambda_star*V,'fro')/ ...
    (max(1,abs(lambda_star))*norm(V,'fro'));
assert(residual<1e-7,'Critical eigenspace residual is too large: %.3g',residual);
mass=sum(abs(V).^2,2)/size(V,2);
weight=sum(reshape(mass,N,K),2); weight=weight/sum(weight);
degree=zeros(N,1);
for k=1:K, degree=degree+full(diag(layers{k}))/K; end
[ranked_degree,order]=sort(degree,'descend');
% Only place vertices at degree-group boundaries, avoiding arbitrary tie order.
boundaries=[find(abs(diff(ranked_degree))>1e-10); N];
cumulative=cumsum(weight(order));
profile=struct('f',[0;boundaries/N], ...
    'cumulative',[0;cumulative(boundaries)],'node_weight',weight, ...
    'mean_node_degree',degree,'lambda_star',lambda_star, ...
    'multiplicity',size(V,2),'residual',residual, ...
    'top10_weight',interp1([0;boundaries/N],[0;cumulative(boundaries)],0.1));
end
