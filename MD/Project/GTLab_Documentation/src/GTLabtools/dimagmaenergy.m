%% Last edited on 25Apr2025

edited='25Apr2025'
%edited='23Apr2025'
%edited='16Apr2025'
%edited='9Jan2025'
switch edited
case '25Apr2025'
% auxiliar functios
% 25Abr2025------------------------------------------------------
% further auxiliar
% activations
        gelu = @(x) 0.5 .* x .* (1 + tanh(sqrt(2/pi) .* (x + 0.044715 .* x.^3)));
        gelu_derivative = @(x) ...
            0.5 .* (1 + tanh(sqrt(2/pi)*(x + 0.044715.*x.^3))) + ...
            0.5 .* x .* (1 - tanh(sqrt(2/pi)*(x + 0.044715.*x.^3)).^2) .* ...
            sqrt(2/pi) .* (1 + 3*0.044715*x.^2);
        softmax = @(z) exp(z - max(z)) ./ sum(exp(z - max(z)));

% auxiliar functions
    function [y3,y2,y1,z3,z2,z1]=forward(x,W1,W2,W3,b1,b2,b3)
        % Forward pass
        nlay=size(W2,3);
        nlin=size(W2,1);
        y2=zeros(nlin,nlay);
        z2=y2;
        z1 = W1 * x + b1;
        y1 = gelu(z1);

        z2(:,1) = W2(:,:,1) * y1 + b2(:,1);
        y2(:,1) = gelu(z2(:,1));
        for i=2:nlay
            z2(:,i) = W2(:,:,i) * y2(:,i-1) + b2(:,i);
            y2(:,i) = gelu(z2(:,i));
        end

        z3 = W3 * y2(:,end) + b3;
        y3 = softmax(z3);  % nB × 1
    end
% function backprop
function [W1, W2, W3, b1, b2, b3] = backprop(x, y_true, y3, y2, y1, z3, z2, z1, W1, W2, W3, b1, b2, b3, eta)
    nlay = size(W2, 3);

    % Output layer
    dz3 = y3 - y_true; % softmax + crossentropy derivative
    dW3 = dz3 * y2(:, end)';  % last y2 layer
    db3 = dz3;

   % Hidden layer 2 backpropagation
dy2 = W3' * dz3;
dz2 = zeros(size(z2));
dW2 = zeros(size(W2));
db2 = zeros(size(b2));

for i = nlay:-1:2
    dz2(:, i) = dy2 .* gelu_derivative(z2(:, i));
    input_to_this_layer = y2(:, i - 1);
    dW2(:, :, i) = dz2(:, i) * input_to_this_layer';
    db2(:, i) = dz2(:, i);
    dy2 = W2(:, :, i)' * dz2(:, i);
end

% Explicit case for i = 1
dz2(:, 1) = dy2 .* gelu_derivative(z2(:, 1));
dW2(:, :, 1) = dz2(:, 1) * y1';
db2(:, 1) = dz2(:, 1);

    % Hidden layer 1
    dy1 = W2(:, :, 1)' * dz2(:, 1);
    dz1 = dy1 .* gelu_derivative(z1);
    dW1 = dz1 * x';
    db1 = dz1;

    % Weight updates
    W1 = W1 - eta * dW1;
    b1 = b1 - eta * db1;

    for i = 1:nlay
        W2(:, :, i) = W2(:, :, i) - eta * dW2(:, :, i);
        b2(:, i) = b2(:, i) - eta * db2(:, i);
    end

    W3 = W3 - eta * dW3;
    b3 = b3 - eta * db3;
end

function Results=show_performance(a_list, a0_list, x0_list, nB, epoch)
    % Shows performance metrics after an epoch:
    % - Number of correct predictions between a and a0
    % - Geometric error ||x0 - x_pred||

    N = size(a_list, 1);  % number of sequences
    nX = size(a_list, 2); % sequence length

    % Define the embedding function (should match your training logic)
    %rho = @(a) 3.^(repmat(1:nX, size(a,1),1) - 1 - floor(nX/2) - floor((a-1)/nB)*nX);
    %geom = @(a) rho(a) .* exp(2i*pi*a/nB);

    % Prediction vs. target symbol match
    correct_symbols = sum(a_list == a0_list, 2);  % correct symbols per sequence
    total_correct = sum(correct_symbols);
    
    % Geometric difference
    x0_pred = sum(geom(a0_list), 2);
    geom_error = abs(x0_list - x0_pred);  % vector of distances

    % Display stats
    fprintf('Epoch %d Results:\n', epoch);
    fprintf('  Total sequences:          %d\n', N);
    fprintf('  Total correct symbols:    %d out of %d (%.2f%%)\n', ...
        total_correct, N * nX, 100 * total_correct / (N * nX));
    fprintf('  Mean geometric error:     %.5f\n', mean(geom_error));
    fprintf('  Max geometric error:      %.5f\n\n', max(geom_error));

    Results=[total_correct / (N * nX) mean(geom_error) max(geom_error)];
end

% the main code starts here
        %case 4 % 25Abr2025
        % a completely different approach
        nX = 3; nB = 4; nA = (3*nB)^nX; % in this interpretation A=(B+B+B)^X

        % i2s and s2i are bijections between A and (B+B+B)^X, so that, in
        % practie, instead of nB we need to use 3*nB
        nB=3*nB;
        i2s = @(i) mod(floor(bsxfun(@rdivide, i(:)-1, nB.^((1:nX)-1))), nB) + 1;
        s2i = @(s) (s-1)*(nB.^((1:nX)'-1)) + 1;
        nB=nB/3; % little trick to avoid changing the definition of s2i and i2s

        % anonymous
        id=@(x) (1:numel(x))';
        % given a=i2s(k), define rho(a) as
        rho=@(a) 3.^(repmat(1:nX,size(a,1),1)-1-floor(nX/2)-floor((a-1)/nB)*nX);

        % to visualize the geometric realization
        a=i2s(1:nA);
        r=rho(a); t=exp(2i*pi*a/nB);
        plot(sum(r.*t,2),'.')

        % embedding
        geom=@(a) rho(a).*exp(2i*pi*a/nB);
        embed=@(x,a) log([x, sum(geom(a),2), geom(a)]);
        split=@(x) [real(x(:)); imag(x(:))];

        % activations
        gelu = @(x) 0.5 .* x .* (1 + tanh(sqrt(2/pi) .* (x + 0.044715 .* x.^3)));
        gelu_derivative = @(x) ...
            0.5 .* (1 + tanh(sqrt(2/pi)*(x + 0.044715.*x.^3))) + ...
            0.5 .* x .* (1 - tanh(sqrt(2/pi)*(x + 0.044715.*x.^3)).^2) .* ...
            sqrt(2/pi) .* (1 + 3*0.044715*x.^2);
        softmax = @(z) exp(z - max(z)) ./ sum(exp(z - max(z)));

        % Network dimensions
        Nepoch=25;
        eta_base = linspace(1e-1,1e-3,Nepoch);
        threshold_macro=.9; threshold_micro=.45; threshold_nano=0.1;
        nW=8; % scale factor for dimensions of W1, W2,W3
        nE= numel(split(embed(1,i2s(nA))));  % dim of embedding
        nW1col =nE;  nW1lin=nW*nE;
        nW2lay=nX; % W2 has nX layers
        nW2col=nW*nE; nW2lin=nW1lin; % W2 is nW*nE-by-nW*nE-by-nX
        nW3col=nW1lin; nW3lin=nB;

        % Initialize weights and biases
        W1 = 1/nE * ones(nW1lin, nW1col);
        b1 = zeros(nW1lin, 1);
        W2 = 1/nW2col * ones(nW2lin, nW2col,nW2lay);
        b2 = zeros(nW2lin, nW2lay);
        W3 = 1/nW3col * ones(nW3lin, nW3col);
        b3 = zeros(nW3lin, 1);
        

        Results=zeros(Nepoch,3);
        for epoch=1:Nepoch
        % training per epoch
        eta=eta_base(epoch);
        a_true=zeros(nA-1,nX); a_pred=zeros(nA-1,nX); x0_list=zeros(nA-1,1);
        for ix=1:nA
            a=i2s(ix); %disp(ix)
            x0=sum(geom(a));
            a0=ones(1,nX); % initial guess for a            
            % first round nano
            for j=1:nX
            x=split(embed(x0,a0));
            [y3,y2,y1,z3,z2,z1]=forward(x,W1,W2,W3,b1,b2,b3);
            [u,v]=max(y3);
            if threshold_nano*(epoch/Nepoch)<u, a0(j)=v+2*nB; end
            y_true=full(sparse(mod(a(j)-1,nB)+1,1,1,nB,1));
            [W1, W2, W3, b1, b2, b3] = backprop(x, y_true, y3, y2, y1, z3, z2, z1, W1, W2, W3, b1, b2, b3, eta);
            end
            % second round micro
            for j=1:nX
            x=split(embed(x0,a0));
            [y3,y2,y1,z3,z2,z1]=forward(x,W1,W2,W3,b1,b2,b3);
            [u,v]=max(y3);
            if threshold_micro*(epoch/Nepoch)<u, a0(j)=v+nB; end
            y_true=full(sparse(mod(a(j)-1,nB)+1,1,1,nB,1));
            [W1, W2, W3, b1, b2, b3] = backprop(x, y_true, y3, y2, y1, z3, z2, z1, W1, W2, W3, b1, b2, b3, eta);
            end
            % third round macro
            for j=1:nX
            x=split(embed(x0,a0));
            [y3,y2,y1,z3,z2,z1]=forward(x,W1,W2,W3,b1,b2,b3);
            [u,v]=max(y3);
            if threshold_macro*(epoch/Nepoch)<u, a0(j)=v; end
            y_true=full(sparse(mod(a(j)-1,nB)+1,1,1,nB,1));
            [W1, W2, W3, b1, b2, b3] = backprop(x, y_true, y3, y2, y1, z3, z2, z1, W1, W2, W3, b1, b2, b3, eta);
            end
            a_true(ix,:)=a; a_pred(ix,:)=a0; x0_list(ix)=x0; %disp(ix)
        end
        Results(epoch,:)=show_performance(a_true, a_pred, x0_list, nB, epoch);
   
        end % epoch
subplot(2,1,1), plot(Results(:,1))
subplot(2,1,2), plot(Results(:,2:3)) 



%------------------------------------------------------------------
case '23Apr2025'
% setting the stage
nX=7; nB=4; nA=nB^nX;

% defining a bijection between A=1:nX^nX and End(X)=X^X
i2s=@(i) mod(floor(bsxfun(@rdivide,i(:)-1,nB.^((1:nX)-1))),nB)+1;
s2i=@(s) (s-1)*(nB.^((1:nX)'-1))+1;

% the auxiliar functions f:X-->Integers and g:B-->Complex
f3=3.^(-ceil((nX-1)/2):floor((nX-1)/2))';
% the definition of f has been changed, instead of f=[0,1,-1,-2,2,-3] we now consider f=[-3,-2,-1,0,1,2]
%f=zeros(nX,1);
%f(2:2:nX)=-(1:ceil((nX-1)/2));
%f(3:2:nX)=(1:floor((nX-1)/2));
g=[exp(2*pi*1i*(nB-(1:nB-1))/(nB-1)), 0];

% defining the map z:A-->Complex, with A=B^X
z=@(a) g(a)*f3;

% our goal is to find the inverse of z, zinv:Complex-->A

% Embedding function

embedding = @(x) [real(x); imag(x); abs(x); angle(x); angle(log(x)); angle(exp(x)); real(sqrt(x)); imag(sqrt(x)); 1+abs(x).^2; 2*real(x)./(1+abs(x).^2); 2*imag(x)./(1+abs(x).^2); (-1+abs(x).^2)./(1+abs(x).^2)];


nE = 11;  % Embedding dimension (updated to match actual feature count)

% Hyperparameters
eta = 1e-6;
nW1 = 20; nW2=30;
nE=numel(embedding(0));
nepoch=1; % number of epochs

% Initialize weights and biases
W1 = 1/nE*ones(nW1, nE);
b1 = zeros(nW1, nB);
W2 = 1/nW1*ones(nW2, nW1);
b2 = zeros(nW2, nB);
W3 = 1/nW2*ones(nX, nW2);
b3 = zeros(nX, nB);


% Activation functions
gelu = @(x) 0.5 .* x .* (1 + tanh(sqrt(2/pi) .* (x + 0.044715 .* x.^3)));
gelu_derivative = @(x) ...
    0.5 .* (1 + tanh(sqrt(2/pi)*(x + 0.044715.*x.^3))) + ...
    0.5 .* x .* (1 - tanh(sqrt(2/pi)*(x + 0.044715.*x.^3)).^2) .* ...
    sqrt(2/pi) .* (1 + 3*0.044715*x.^2);

%softmax = @(z) exp(z) ./ sum(exp(z));
softmax = @(z) exp(z - max(z)) ./ sum(exp(z - max(z))); % more stable
cross_entropy_loss = @(y_pred, y_true) -1/nX*sum(sum(y_true .* log(y_pred + 1e-12)));  % Stability

% Tracking
loss_vals = zeros(nA, 1);
accuracy_vals = zeros(nA, 1);

% Initialize a counter for correct predictions
correct_preds = 0;

% Training loop
for i=1:nepoch
for k = 1:nA
    y = i2s(k); % One-hot vector (1-by-nX)
    x0 = z(y) + f3(1)/2*g;
 y_true=full(sparse(1:nX,y,1,nX,nB));
    
    % Apply embedding transformation
    x = embedding(x0);  % Use the further embedding function here

    % Forward pass
    z1 = W1 * x + b1;
    y1 = gelu(z1);

    z2 = W2 * y1 + b2;
    y2 = gelu(z2);

    z3 = W3 * y2 + b3;
    y3 = softmax(z3);  % Final layer: Probabilities for each class (7x7)
    if any(isnan(y3)), disp(k), break, end
 % Loss
    loss_vals(k) = cross_entropy_loss(y3, y_true);  % Use full vector for cross-entropy

    % Accuracy calculation
    [u, y_pred] = max(y3, [], 2);  % Predicted class (index of max probability)
    correct_preds =  nnz(y_pred == y');  % Accumulate correct predictions
    
    % Update accuracy at each iteration
    accuracy_vals(k) = correct_preds/nX;  

    % Backpropagation
    dz3 = y3 - y_true;
    dW3 = dz3 * y2';
    db3 = dz3;

    dy2 = W3' * dz3;
    dz2 = dy2 .* gelu_derivative(z2);
    dW2 = dz2 * y1';
    db2 = dz2;

    dy1 = W2' * dz2;
    dz1 = dy1 .* gelu_derivative(z1);
    dW1 = dz1 * x';
    %db1 = db1 + dz1;
    db1 =  dz1;
    
    % Update weights
    W1 = W1 - eta * dW1;
    b1 = b1 - eta * db1;
    W2 = W2 - eta * dW2;
    b2 = b2 - eta * db2;
    W3 = W3 - eta * dW3;
    b3 = b3 - eta * db3;

   % Display progress every 100 samples
    if mod(k, 1000) == 0 || k == 1
        fprintf('Iter %d | Loss: %.4f | Accuracy: %d\n', ...
                k, loss_vals(k), correct_preds);
                disp([y_pred y'])
%disp([ y3 y_true])
    end

end
end % nepoch

% Plotting
figure;
subplot(2,1,1); plot(cumsum(loss_vals/nA)); title('Cross-Entropy Loss'); xlabel('Sample'); ylabel('Loss');
subplot(2,1,2); plot(cumsum(accuracy_vals/nX)); title('Accuracy per Sample'); xlabel('Sample'); ylabel('Accuracy');

% end of case 23Apr2025



%------------------------------------
case '16Apr2025'
% setting the stage
nX=6; nB=nX; nA=nB^nX;

% defining a bijection between A=1:nX^nX and End(X)=X^X
i2s=@(i) mod(floor(bsxfun(@rdivide,i(:)-1,nB.^((1:nX)-1))),nB)+1;
s2i=@(s) (s-1)*(nB.^((1:nX)'-1))+1;

% the auxiliar functions f:X-->Integers and g:B-->Complex
f=zeros(nX,1);
f(2:2:nX)=-(1:ceil((nX-1)/2));
f(3:2:nX)=(1:floor((nX-1)/2));
g=[exp(2*pi*1i*(nB-(1:nB-1))/(nB-1)), 0];

% defining the map z:A-->Complex, with A=B^X
z=@(a) g(a)*(3.^f);

% our goal is to find the inverse of z, zinv:Complex-->A

% Embedding function
%complex2sphere = @(x) [real(x)/abs(x); imag(x)/abs(x)]; % Define or load if needed
embedding = @(x) [real(x); imag(x); abs(x); angle(x); abs(log(x)); angle(log(x)); ...
                  real(exp(x)); imag(exp(x)); complex2sphere(x)'];
furtherembbed=@(x) embedding(x)*1./(1:nX);

nE = 11;  % Embedding dimension (updated to match actual feature count)

% Hyperparameters
eta = 1e-6;
nW = 20;

% Initialize weights and biases
W1 = randn(nW, nE);
b1 = randn(nW, nX);
W2 = randn(nW, nW);
b2 = randn(nW, nX);
W3 = randn(nX, nW);
b3 = randn(nX, nX);


% Activation functions
gelu = @(x) 0.5 .* x .* (1 + tanh(sqrt(2/pi) .* (x + 0.044715 .* x.^3)));
gelu_derivative = @(x) ...
    0.5 .* (1 + tanh(sqrt(2/pi)*(x + 0.044715.*x.^3))) + ...
    0.5 .* x .* (1 - tanh(sqrt(2/pi)*(x + 0.044715.*x.^3)).^2) .* ...
    sqrt(2/pi) .* (1 + 3*0.044715*x.^2);

%softmax = @(z) exp(z) ./ sum(exp(z));
softmax = @(z) exp(z - max(z)) ./ sum(exp(z - max(z))); % more stable
cross_entropy_loss = @(y_pred, y_true) -1/nX*sum(sum(y_true .* log(y_pred + 1e-12)));  % Stability

% Tracking
loss_vals = zeros(nA, 1);
accuracy_vals = zeros(nA, 1);

% Initialize a counter for correct predictions
correct_preds = 0;

% Training loop
for k = 1:nA
    y = i2s(k); % One-hot vector (7x1)
    x0 = z(y) + 0.5 * 3^(-nX/2) * exp(randi(360) * 2 * pi * 1i);
 y_true=full(sparse(1:nX,y,1,nX,nX));
    
    % Apply further embedding transformation
    x = furtherembbed(x0);  % Use the further embedding function here

    % Forward pass
    z1 = W1 * x + b1;
    y1 = gelu(z1);

    z2 = W2 * y1 + b2;
    y2 = gelu(z2);

    z3 = W3 * y2 + b3;
    y3 = softmax(z3);  % Final layer: Probabilities for each class (7x7)
    if any(isnan(y3)), break, end
 % Loss
    loss_vals(k) = cross_entropy_loss(y3, y_true);  % Use full vector for cross-entropy

    % Accuracy calculation
    [~, pred_class] = max(y3, [], 2);  % Predicted class (index of max probability)
    correct_preds =  sum(pred_class == y');  % Accumulate correct predictions
    
    % Update accuracy at each iteration
    accuracy_vals(k) = correct_preds / nX;  % Divide by the total samples processed so far

    % Backpropagation
    dz3 = y3 - y_true;
    dW3 = dz3 * y2';
    db3 = dz3;

    dy2 = W3' * dz3;
    dz2 = dy2 .* gelu_derivative(z2);
    dW2 = dz2 * y1';
    db2 = dz2;

    dy1 = W2' * dz2;
    dz1 = dy1 .* gelu_derivative(z1);
    dW1 = dz1 * x';
    %db1 = db1 + dz1;
    db1 =  dz1;
    
    % Update weights
    W1 = W1 - eta * dW1;
    b1 = b1 - eta * db1;
    W2 = W2 - eta * dW2;
    b2 = b2 - eta * db2;
    W3 = W3 - eta * dW3;
    b3 = b3 - eta * db3;

   % Display progress every 100 samples
    if mod(k, 1000) == 0 || k == 1
        fprintf('Iter %d | Loss: %.4f | Accuracy: %.2f%%\n', ...
                k, loss_vals(k), accuracy_vals(k) * 100);
                disp([pred_class y'])
%disp([ y3 y_true])
    end

end


% Plotting
figure;
subplot(2,1,1); plot(loss_vals); title('Cross-Entropy Loss'); xlabel('Sample'); ylabel('Loss');
subplot(2,1,2); plot(accuracy_vals); title('Accuracy per Sample'); xlabel('Sample'); ylabel('Accuracy');

% end of case 16Apr2025

%----------------------------------
case '9Jan2025'
disp("On the energy of dimagmas - Math4Games 9Jan2025");
% the case N=2, M=3
nX=2 % N=2
nB=3 % M=3
X=1:nX
B=1:nB
A=fillzeros(zeros(nX,1),B) % use function fillzeros.m available on Teams->Files->Matlab-Octave 
nA=size(A,2)
% check that nA=nB^nX
isequal(nA,nB^nX)
% construct the map f:X-->Integers
f=zeros(1,nX); % inicialize f with zeros
f(2:2:end)=-(2:2:nX)/2; % f(x) if x is even
f(1:2:end)=((1:2:nX)-1)/2 % f(x) if x is odd
% construct the map g:B-->Complex
g=[exp(2*pi*1i*(nB-B(1:end-1))/(nB-1)), 0]
% construct the map F:A-->Complex
F=@(alpha) sum(g(alpha).*3.^repmat(f',1,size(alpha,2))) % use anonimous function 
% compute F(A), recall that each column of A encodes one particular alpha
Fall=F(A) % Fall contains all values of F(alpha), for all possible alphas
% plot Fall
figure(1), clf
plot(Fall,'--'), ax=max(abs(Fall))+1; axis([-ax ax -ax ax])
text(real(Fall),imag(Fall),num2str((1:nA)'))
% check if F is injective, i.e., Fall has no repeated elements
if any(abs(diff(sort(Fall)))<1e-12) % 12 significant digits
disp('Up to 12 digits, F is not injective')
else
    disp('F is injective')
end
% generating Plus and Times tables using modular arithmetic
[i,j]=ndgrid(1:nA,1:nA)
% note that any other nA-by-nA two tables with values in the range 1:nA would do as well, this is just one example
Plus=mod(i+j,nA)+1 
Times=mod(i.*j,nA)+1 
% Finally, The Energy of the dimagma (Plus,Times)
DiffPlus=Fall(i(:))+Fall(j(:))-Fall(Plus(sub2ind([nA,nA],i(:),j(:))));
DiffTimes=Fall(i(:)).*Fall(j(:))./Fall(Times(sub2ind([nA,nA],i(:),j(:))))-1;
DiffTimes(~isfinite(DiffTimes))=0;
figure(2)
subplot(1,2,1)
plot(sort(abs(DiffPlus)),'.-')
subplot(1,2,2)
plot(sort(abs(DiffTimes)),'.-')
E=sum(abs(DiffPlus))+sum(abs(DiffTimes))

%------------------------------------------
%%%%%  The Case N=4, M=9 %%%%%%%%%%%%%%%
% Not showing outputs for obvious reasons,
% appart form that, the code is the same as before
% no further comments
% the case N=4, M=9
nX=4 % N=4, 2
nB=9 % M=9, 25
X=1:nX
B=1:nB
A=fillzeros(zeros(nX,1),B);nA=size(A,2);
isequal(nA,nB^nX)
f=zeros(1,nX); % inicialize f with zeros
f(2:2:end)=-(2:2:nX)/2; % f(x) if x is even
f(1:2:end)=((1:2:nX)-1)/2; % f(x) if x is odd
g=[exp(2*pi*1i*(nB-B(1:end-1))/(nB-1)), 0];
F=@(alpha) sum(g(alpha).*3.^repmat(f',1,size(alpha,2)));
Fall=F(A); % Fall contains all values of F(alpha), for all possible alphas, i.e., all alpha_i with i ranging 1:nA
figure(3)
plot(Fall,'.'), ax=max(abs(Fall))+1; axis([-ax ax -ax ax])
%text(real(Fall),imag(Fall),num2str((1:nA)'))
% check if F is injective, i.e., Fall has no repeated elements
if any(abs(diff(sort(Fall)))<1e-12) % 12 significant digits
disp('Up to 12 digits, F is not injective')
else
    disp('F is injective')
end
% generating Plus and Times tables using modular arithmetic
if nA^2<1e6 % otherwise we run into problems
[i,j]=ndgrid(1:nA,1:nA);
% note that any other nA-by-nA two tables with values in the range 1:nA would do as well, this is just one example
Plus=mod(i+j,nA)+1; 
Times=mod(i.*j,nA)+1; 
% Finally, The Energy of the dimagma (Plus,Times)
DiffPlus=Fall(i(:))+Fall(j(:))-Fall(Plus(sub2ind([nA,nA],i(:),j(:))));
DiffTimes=Fall(i(:)).*Fall(j(:))./Fall(Times(sub2ind([nA,nA],i(:),j(:))))-1;
DiffTimes(~isfinite(DiffTimes))=0;
E=sum(abs(DiffPlus))+sum(abs(DiffTimes))*1i
else
    disp('Too big to compute the energy in one stroke')
end

% To run several cases in one shot
clear, clf, cla
%---------------------------------------
k=0; N=2:4; M=3:6;
for k_N=1:numel(N)
for k_M=1:numel(M)
k=k+1
nX=N(k_N); nB=M(k_M);
X=1:nX;
B=1:nB;
A=fillzeros(zeros(nX,1),B);nA=size(A,2);
isequal(nA,nB^nX)
f=zeros(1,nX); % inicialize f with zeros
f(2:2:end)=-(2:2:nX)/2; % f(x) if x is even
f(1:2:end)=((1:2:nX)-1)/2; % f(x) if x is odd
g=[exp(2*pi*1i*(nB-B(1:end-1))/(nB-1)), 0];
F=@(alpha) sum(g(alpha).*3.^repmat(f',1,size(alpha,2)));
Fall=F(A); % Fall contains all values of F(alpha), for all possible alphas, i.e., all alpha_i with i ranging 1:nA
figure(1)
subplot(4,3,k)
plot(Fall,'.'),  axis off
title(['F; N=', num2str(nX) ', M= ', num2str(nB)])
figure(2)
subplot(4,3,k)
plot(log(Fall),'.'),  axis off
title(['log(F); N=', num2str(nX) ', M= ', num2str(nB)])
figure(3)
subplot(4,3,k)
plot(exp(Fall),'.'),  axis off
title(['exp(F); N=', num2str(nX) ', M= ', num2str(nB)])
end
end

end % edited switch case
