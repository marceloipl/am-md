disp("Aula MDPL25 dia 7Mar2025");

%class=1; %7Mar
class=2; %20Mar
switch class
case 1 % class 7Mar 
% a first implementation of fillzeros
A0=[0 2 0 -1]
A1=A0(:)
z=find(A1==0)  % z=[1 3]'
nz=numel(z)
[a b]=ndgrid([1 3], 1:3)
ab=[a(:) b(:)]
A=repmat(A1,1,size(ab,1))
A(z,:)=ab'
%{
for i=1:nz
    A(A1==-i,:)=A(z(i),:);
end
% this has a problem when i=2 because on the left hand side of A(A1==-i,:)=A(z(i),:) we get an empty matrix whereas on the right hand side we get the the third line of A; a better solution is:
%}
A(A1<0,:)=A(z(abs(A1(A1<0))),:);
% replacing negative values by the respective counterparts

% a first study on determining the number of conjugacy  classes of endomaps on a three-element set 

% A is the collection of all endomaps of X=1:3
A=fillzeros([0 0 0]',1:3)
Circ=zeros(27,27);
for i=1:27
for j=1:27
ui=A(:,i);
uj=A(:,j);
uij=ui(uj);
Circ(i,j)=find(ismember(A',uij(:)','rows'));
end
end

% identify the neutral element:
% (a) in the index matrix A; (b) in the multiplication table Circ
find(ismember(A',[1 2 3],'rows'))
find(ismember(Circ,1:27,'rows'))

% identify the invertible elements (i.e. permutations):
% (a) in the index matrix A; (b) in the multiplication table Circ
ismember(sort(A)',[1 2 3],'rows')'
ismember(sort(Circ,2),1:27,'rows')'
% disp(A)

% another possibility is
Perm=perms(1:3)'
Inv=find(ismember(A',Perm','rows'))

% Define the relation of equivalence up to conjugation, i.e., for all u and v in A, we say that u is related to v (and write uRv) iff there exists p in Per such that p(u)=v(p)

R=zeros(27,27);
for i=1:27
for j=1:27
for k=1:6
u=A(:,i);
v=A(:,j);
p=Perm(:,k);
if p(u)==v(p), R(i,j)=1; end
end
end
end

% to extract the graph associated with the relation R we do:
[d,c]=find(R);

% And the quotient A/R is computed with the coequalizer
q=coeq(d,c)

% The number of equivalence classes is:
numel(unique(q))


case 2 % class 20Mar
%%%%-----------------------------------
% Energy of dimagmas with neural network optimization

clear, clf, cla
nX=3; nB=nX; nA=nB^nX; X=1:nX; B=1:nB; % in this case B isequal X
A=(1:nA)';

% defining a bijection between A=1:nX^nX and End(X)=X^X

i2s=@(i) mod(floor(bsxfun(@rdivide,i(:)-1,nB.^((1:nX)-1))),nB)+1;
s2i=@(s) (s-1)*(nB.^((1:nX)'-1))+1;

% defining the maps f and g
f=zeros(1,nX); % inicialize f with zeros
f(2:2:end)=-(2:2:nX)/2; % f(x) if x is even
f(1:2:end)=((1:2:nX)-1)/2; % f(x) if x is odd
g=[exp(2*pi*1i*(nB-B(1:end-1))/(nB-1)), 0];

R=sum(3.^f);
h=sqrt(2)/(2*nA);
zaux=@(a) g(a).*(3.^repmat(f,size(a,1),1));
z=@(a) sum(zaux(a),2);
w=@(a) exp(-imag(z(a))+1i*pi*(real(z(a))/R));
F=[z(i2s(A)); w(i2s(A))];

plot(F,'.')

% Using neural networks to find an addition at the level of the indexes, using only z and zaux, leaving w for a later study

% defining an activation map tau:Reals-->Reals
M=nB;  % M=3
tau=@(t) (1+(M-2)/2*(cos(t)+1)).*(t>0) + M*(t<=0);
tauprime=@(t) -(M-2)/2*sin(t).*(t>0);
plot(linspace(-10,10),tau(linspace(-10,10)),'.-')
axis([-10 10 0 M])

% initializing parameters
W=ones(3,9); b=ones(3,1); % it could be random

Energy=0;
for k=1:200
a1=i2s(randi(nA)); a2=i2s(randi(nA)); a3=i2s(randi(nA));
y=z(a1)+z(a2)*(z(a3)-z(a1)); % y is the exact value
x=[g(a1) g(a2) g(a3)]'; % x is the embedding of a1 and a2

% forward
yhat=round(tau(real(W*x+b)))'; % yhat is the computed value to be compared with y in the loss function

%disp([a1 i2s(a1) a2 i2s(a2) yhat])
L=(z(yhat)-y).^2; % the loss function
Energy=Energy+sqrt(abs(L));
%disp([sqrt(abs(L)) Energy])
disp([real([z(a1) z(a2) z(a3)]) real(y) real(z(yhat))])  % it is supposed to get smaller form one iteration to another
% backward
eta=0.1;  % a parameter that can be adjusted
partialder=[(2*pi*1i*(nB-B(1:end-1))/(nB-1)) 0]'.*tauprime(real(W*x+b));
b=b-eta*2*(z(yhat)-y)*zaux(yhat)'.*partialder;
W=W-eta*2*(z(yhat)-y)*(x*zaux(yhat))'.*repmat(partialder,1,9);
end

% testing the total accuracy of the model
Energy=0;
for i1=1:nA
    for i2=1:nA
        for i3=1:nA
            a1=i2s(i1); a2=i2s(i2); a3=i2s(i3);
        y=z(a1)+z(a2)*(z(a3)-z(a1)); % y is the exact value
        yhat=round(tau(real(W*x+b)))';
        L=(z(yhat)-y).^2;
        Energy=Energy+sqrt(abs(L));
    end
end
end
disp(Energy)

%---------------------------------------------
%% Another trial, as suggested by ChatGPT 
% --------------------------------------------
% A generic implementation of a neural network with 2 layers and sigmoid activation function

% Sigmoid activation function and its derivative
sigmoid = @(x) 1 ./ (1 + exp(-x));
sigmoid_derivative = @(x) sigmoid(x) .* (1 - sigmoid(x));

% Loss function (Mean Squared Error)
mse_loss = @(y_pred, y_true) 0.5 * sum((y_pred - y_true).^2);

% Forward pass function
function [z1, y1, z2, y2] = forward_pass(x, W, b, W_prime, b_prime, sigmoid)
    z1 = W * x + b;
    y1 = sigmoid(z1);
    z2 = W_prime * y1 + b_prime;
    y2 = sigmoid(z2);
end

% Backward pass function
function [W, b, W_prime, b_prime] = backward_pass(x, y_true, z1, y1, z2, y2, W, b, W_prime, b_prime, eta,sigmoid_derivative)
    % Compute loss gradient with respect to y2
    dL_dy2 = y2 - y_true;
    
    % Gradients for second layer (W' and b')
    dL_dz2 = dL_dy2 .* sigmoid_derivative(z2);
    dL_dW_prime = dL_dz2 * y1';
    dL_db_prime = dL_dz2;
    
    % Gradients for first layer (W and b)
    dL_dy1 = W_prime' * dL_dz2;
    dL_dz1 = dL_dy1 .* sigmoid_derivative(z1);
    dL_dW = dL_dz1 * x';
    dL_db = dL_dz1;
    
    % Update parameters using gradient descent
    W_prime = W_prime - eta * dL_dW_prime;
    b_prime = b_prime - eta * dL_db_prime;
    W = W - eta * dL_dW;
    b = b - eta * dL_db;
end

% testing with a trivial example

% Example input and true output
x = [0.5; -0.5];
y_true = [0.1];

% Initialize weights and biases
W = randn(1, 2);
b = randn(1, 1);
W_prime = randn(1, 1);
b_prime = randn(1, 1);

% Learning rate
eta = 0.01;

% Perform forward pass
[z1, y1, z2, y2] = forward_pass(x, W, b, W_prime, b_prime,sigmoid);

% Perform backward pass and update parameters
[W, b, W_prime, b_prime] = backward_pass(x, y_true, z1, y1, z2, y2, W, b, W_prime, b_prime, eta,sigmoid_derivative);

% Print the updated weights and biases
disp('Updated W:');
disp(W);
disp('Updated b:');
disp(b);
disp('Updated W_prime:');
disp(W_prime);
disp('Updated b_prime:');
disp(b_prime);

% testing with our example of dimagma energy on the complex plane with nX=nB=3
nX=3; nB=nX; nA=nB^nX; X=1:nX; B=1:nB; % in this case B isequal X
A=(1:nA)';

% defining a bijection between A=1:nX^nX and End(X)=X^X
i2s=@(i) mod(floor(bsxfun(@rdivide,i(:)-1,nB.^((1:nX)-1))),nB)+1;
s2i=@(s) (s-1)*(nB.^((1:nX)'-1))+1;

% defining the maps f and g
f=zeros(1,nX); % inicialize f with zeros
f(2:2:end)=-(2:2:nX)/2; % f(x) if x is even
f(1:2:end)=((1:2:nX)-1)/2; % f(x) if x is odd
g=[exp(2*pi*1i*(nB-B(1:end-1))/(nB-1)), 0];

zaux=@(a) g(a).*(3.^repmat(f,size(a,1),1));
z=@(i) sum(zaux(i2s(i)),2);

% to extract the index that best approximates a given value of z 
[zord,iord]=sort(real(z(1:nA)));
%iordinv(iord)=1:nA;
tauinv=@(z) iord(round(interp1(zord,1:nA,real(z))));
% tauinv(z(i)) is i
isequal(tauinv(z(1:nA)),(1:nA)')

[i1 i2 i3]=ndgrid(1:nA);
i123=[i1(:) i2(:) i3(:)];
z123=z(i1)+z(i2).*(z(i3)-z(i2));
I_=[i123, real(z123(:))];
I=I_(abs(z123(:))<4.5,1:3);

% the iterative process

% Initialize weights and biases
W = randn(3, 9);
b = randn(3, 1);
W_prime = randn(3, 9);
b_prime = randn(3, 1);

% Learning rate
eta = 0.01;

% to be sorted out
for k=1:200
a_=I(randi(size(I,1)),:);
a1=i2s(a_(1)); a2=i2s(a_(2)); a3=i2s(a_(3));
y=z(a_(1))+z(a_(2))*(z(a_(3))-z(a_(1))); % y is the exact value
y_true=g(i2s(tauinv(y)));
x=[g(a1) g(a2) g(a3)]'; % x is the embedding of a1 and a2

% Perform forward pass
[z1, y1, z2, y2] = forward_pass(x, W, b, W_prime, b_prime,sigmoid);

% Perform backward pass and update parameters
[W, b, W_prime, b_prime] = backward_pass(x, y_true, z1, y1, z2, y2, W, b, W_prime, b_prime, eta,sigmoid_derivative);

end

% an alternative, simpler definition of z:A-->Complex and its inverse
f=zeros(nX,1);
f(2:2:nX)=-(1:ceil((nX-1)/2));
f(3:2:nX)=(1:floor((nX-1)/2));
g=[exp(2*pi*1i*(nB-(1:nB-1))/(nB-1)), 0];

z=@(i) g(i2s(i))*(3.^f);
% the inverse of z, zinv:Complex-->A
nA=nB^nX;
% given z0 in Complex, we obtain iz0 such that z(iz0) is the nearest point to z0 indexed by 1:nA
z0=2.5
[~,iz0]=min(abs(z(1:nA)-z0));

end % switch
