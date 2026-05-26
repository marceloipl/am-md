function A=fillzeros(A0,varargin)
% A=fillzeros(A0,varargin)
% 15Fev2025
% Given a matrix A0 with nz zero entries and possibly negative numbers -1,-2,...,-nz in several positions, returns a matrix A, with one more dimension than A0, along which each zero element in A0 is expanded by the input  values in varargin. Each one of the negative entries is then replaced by the values generated in each one of the respective positions of the respective initial zero entries.
% Example:
% A=fillzeros([1 -1 -2; 0 1 -3; 0 0 1 ],2:3,4:5,6:7)
% Further examples at eof.

A1=A0(:); %linearize the matirx A0
z=find(A1==0); % find zeros
nz=numel(z); % count the number of zeros and control the input varargin
assert(length(varargin) == 1 || length(varargin) == nz);
%a=ndgridarr(nz,values);
grid_cells = cell(1, nz);
[grid_cells{:}] = ndgrid(varargin{:});
cab = reshape(cat(3,grid_cells{:}),[],nz);
A=repmat(A1,[1,size(cab,1)]);
A(z,:)=cab';
% replacing negative values by the respective counterparts
A(A1<0,:)=A(z(abs(A1(A1<0))),:);
% replacing negative values by the respective counterparts
A(A1<0,:)=A(z(abs(A1(A1<0))),:);
% returning to the original size of A0
A=reshape(A,[size(A0),size(cab,1)]);
A=squeeze(A);


%{
 A=fillzeros([0 0 ]',1:3)
   1   2   3   1   2   3   1   2   3
   1   1   1   2   2   2   3   3   3
A=fillzeros([0 0 ]',1:3,1:2)
   1   2   3   1   2   3
   1   1   1   2   2   2
A=fillzeros([1 -1 6; 0 1 -2; 6 0 1 ],2:3,4:5)
ans(:,:,1) =

   1   2   6
   2   1   4
   6   4   1

ans(:,:,2) =

   1   3   6
   3   1   4
   6   4   1

ans(:,:,3) =

   1   2   6
   2   1   5
   6   5   1

ans(:,:,4) =

   1   3   6
   3   1   5
   6   5   1

%}

function grid_array = ndgridarr(n, varargin)
% https://stackoverflow.com/questions/60714394/convert-output-of-ndgrid-to-a-single-array
    
    assert(length(varargin) == 1 || length(varargin) == n);

    grid_cells = cell(1, n);
    [grid_cells{:}] = ndgrid(varargin{:});
    %grid_array = reshape(cat(n+1,grid_cells{:}),[],n);
    %grid_array = cell2mat(cellfun(@(c) c(:), grid_cells, 'UniformOutput', false));
    grid_array = reshape(cat(3,grid_cells{:}),[],n);



%----------------------------------------------------------------
% Older versions:-


%{
% 11Fev2024
% Given a matrix A0 with some zero entries, returns a matrix A with one more dimension, along which each zero element in A0 is expanded by the input  values in varargin
% Examples:
% fillzeros([0 0 ]',1:3)
%   1   2   3   1   2   3   1   2   3
%   1   1   1   2   2   2   3   3   3
%fillzeros([0 0 ]',1:3,1:2)
%   1   2   3   1   2   3
%   1   1   1   2   2   2

A1=A0(:); %linearize the matirx A0
z=(A1<=0);
nz=sum(z); 
assert(length(varargin) == 1 || length(varargin) == nz);
%a=ndgridarr(nz,values);
grid_cells = cell(1, nz);
[grid_cells{:}] = ndgrid(varargin{:});
a = reshape(cat(3,grid_cells{:}),[],nz);
A=repmat(A1,[1,size(a,1)]);
A(z,:)=A(z,:)+a';
% returning to the original zize of A0
A=reshape(A,[size(A0),size(a,1)]);
A=squeeze(A);
%}


%{
% before 11Fev2024
% A=fillzeros(A0,values)
% Dec2024 - now using ndgridarr.m instead of ntuples.m 
% the behaviour is still  similar to old fillzeros but allowing parametrized values
% A0 is a matrix with positive integer values and possible with 0 at some
% positions; returns a matrix A with one more dimension than A0 where the
% values 0 are replaced by all possible combinations of input values along the extra dimension
%
% Example
% A=fillzeros([0 1 0]',3:5)
% 
% A =
% 
%      3     4     5     3     4     5     3     4     5
%      1     1     1     1     1     1     1     1     1
%      3     3     3     4     4     4     5     5     5
% 

A1=A0(:); %linearize the matirx A0
%z=(A1==0); % as before
z=(A1<=0);
nz=sum(z); 
%values=1:max(A1);
a=ntuples(nz,values); % as before
%a=ndgridarr(nz,values); % to be explored
A=repmat(A1,[1,size(a,1)]);
%A(z,:)=a';  % as before
A(z,:)=A(z,:)+a';
% returning to the original zize of A0
A=reshape(A,[size(A0),size(a,1)]);
A=squeeze(A);
%}


%-------------------------------------------------------------------------
function A=ntuples(n,values)
% A=ntuples(n,values)
% returns a matrix A whose rows consist of all the possible n-tuples build
% with entries from the vector values
%
% Example
%>> ntuples(3,[4 5])
%
%ans =
%
%     4     4     4
%     5     4     4
%     4     5     4
%     5     5     4
%     4     4     5
%     5     4     5
%     4     5     5
%     5     5     5
%

% the first step is to build the string
% [a1,a2,...,an]=ndgrid(values)
% and evaluate it

% numel(values)^n should be less than 2e7
if numel(values)^n > (2e7), error(['numel(values)^n = ' num2str(numel(values)^n) ' is greater than 20 milions (2e7)']), end

str='[';
for i=1:n
    str=[str 'a' num2str(i) ','];  % last comma has to be ignired, str(1:end-1)
end
str=[str(1:end-1) ']=ndgrid(values);'];
eval(str); % eval('[a1,a2,...,an]=ndgrid(values);')

% the second step is to eval the string
% A=[a1(:), a2(:), ... an(:)];
str='A=[';
for i=1:n
    str=[str 'a' num2str(i) '(:),'];  % last comma has to be ignired, str(1:end-1)
end
str=[str(1:end-1) '];'];
eval(str); % eval('A=[a1(:), a2(:), ... an(:)];')
