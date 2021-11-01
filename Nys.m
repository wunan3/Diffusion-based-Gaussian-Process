function [A] =Nys(data1,data2,data3,y,k,eps,t,sig)
% data1:labeled points over the domain S;
% data2: unlabeled points over the domain S;
% data3: a large data set over the domain S over which we want to extend
% our prediction by Nys;
% y: response variables over labeled points(data1);
% k ,eps, t, sig: GL-GP covariance parameters;
% A: output, the extension of the prediction over data3.

data=[data1;data2];
datae=[data1;data3];
[n, ~] = size(data) ;
[n1, ~] = size(data1) ;
[n2, ~] = size(data2) ;
[n3, ~] = size(data3) ;
[ne, ~] = size(datae) ;

[index,distance]= knnsearch(data, data,'k', n);
distance(:,1)=0;
ker = exp(-distance.^2/(4*eps));
ii = (1:n)'*ones(1,n);
W = sparse(ii, index, ker, n, n);
D = sum(W, 2);
W = bsxfun(@rdivide, bsxfun(@rdivide, W, D), transpose(D));
D = sqrt(sum(W, 2)); 
W = bsxfun(@rdivide, bsxfun(@rdivide, W, D), transpose(D));
[U,S] = eigs(W, k);
U = bsxfun(@rdivide, U(:, 1:end), D);

ker = exp(-distance.^2/(4*eps));
ii = (1:n)'*ones(1,n);
W = sparse(ii, index, ker, n, n); W=full(W);
D1 = sum(W, 2); D1=1./D1;D1=diag(D1);

[index1,distance1]= knnsearch(data, datae,'k', ne);
ker=exp(-distance1.^2/(4*eps));
ii = (1:ne)'*ones(1,n);
W = sparse(ii, index1, ker); W=full(W);
D2 = sum(W, 2); D2=1./D2; D2=diag(D2);
W=D2*W*D1;
D=sum(W, 2); D=1./D; D=diag(D);
E=D*W;


% Construct heat kernel
C = zeros(n,n);
for i=1:k
C = C + exp(t*(S(i,i)-1)/eps)/S(i,i)^2*U(1:n,i)*U(1:n,i)'/norm(U(1:n,i))^2*n;
end
% Prediction
C=E*C*E';
A=mean(y)*ones(n3,1)+C(n1+1:ne,1:n1)/(C(1:n1,1:n1)+sig*eye(n1))*(y-mean(y)*ones(n1,1));


