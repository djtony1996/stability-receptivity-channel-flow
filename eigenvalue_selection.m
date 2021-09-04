% This code is used to construct the truncated energy growth rate matrix.
 
function nA = eigenvalue_selection(A,IMIN,C)

[V,D]      = eig(A);                
D          = diag(D);
[~,is]     = sort(real(D),'descend');
VS         = V(:,is); 
DS         = D(is);

k          = find(real(DS)<IMIN,1);
DS(k:end)  = [];
VS(:,k:end) = [];

WORK = (C*VS)' * C*VS;

[U,S,V] = svd(WORK);
s       = sqrt(diag(S));
F       = diag(s)*U';
invF    = U*diag(ones(size(s))./s);

nA = F * diag(DS) * invF;


end