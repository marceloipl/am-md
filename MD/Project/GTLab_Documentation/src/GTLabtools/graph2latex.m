function s = graph2latex(d, c, V, varargin)
    % GRAPH2LATEX Generate LaTeX TikZ code from a directed graph
    %   s = graph2latex(d, c, V, vL, eL)
    %   - d, c: vectors with indices of domain and codomain (edges)
    %   - V: Nx1 (complex) or Nx2 (real) matrix of node coordinates
    %   - vL (optional): numeric labels for nodes
    %   - eL (optional): numeric labels for edges

    n = size(V, 1);      % Number of nodes
    m = numel(d);        % Number of edges

    % Handle optional arguments
    vL = 1:n;
    eL = [];

    if numel(varargin) >= 1
        vL = varargin{1};
    end
    if numel(varargin) >= 2
        eL = varargin{2};
    end

    % Determine coordinates
    if isvector(V) && ~isreal(V)
        x = real(V);
        y = imag(V);
    elseif size(V, 2) == 2
        x = V(:, 1);
        y = V(:, 2);
    else
        error('V must be Nx1 complex or Nx2 real matrix of coordinates');
    end

    % Header
    header = {
        '\begin{center}'
        '\begin{tikzpicture}[scale=1.5, every node/.style={circle,draw}, every path/.style={->, thick}]'
    };

    % Nodes
    nodes = cell(n, 1);
    for i = 1:n
        labelStr = num2str(vL(i));
        nodes{i} = sprintf('\\node (n%d) at (%.2f, %.2f) {%s};', i, x(i), y(i), labelStr);
    end

    % Edges
    edges = cell(m, 1);
    for i = 1:m
        if isempty(eL)
            edges{i} = sprintf('\\path (n%d) edge (n%d);', d(i), c(i));
        else
            labelStr = num2str(eL(i));
            edges{i} = sprintf('\\path (n%d) edge node[midway, above] {%s} (n%d);', d(i), labelStr, c(i));
        end
    end

    % Footer
    footer = {
        '\end{tikzpicture}'
        '\end{center}'
    };

    % Concatenate all lines
    lines = [header; nodes; edges; footer];
    s = strjoin(lines, '\n');

    % Output to console
    fprintf('%s\n', s);
end

%{
% Example
V = [1+1i; 2+2i; 3+1i; 2+3i];
d = [1; 2; 3; 4];
c = [2; 3; 4; 1];
vL = [10, 20, 30, 40];
eL = [100, 200, 300, 400];

latex_code = graph2latex(d, c, V, vL, eL);

[d,c,v]=hasse(1008)
div=@(n) find(mod(n,1:n)==0);
% and respective Hasse diagram
u=div(1008); % [x,y]=ndgrid(id(u));R=(mod(u(y),u(x))==0); 

str=graph2latex(d,c,v,u);

%}
