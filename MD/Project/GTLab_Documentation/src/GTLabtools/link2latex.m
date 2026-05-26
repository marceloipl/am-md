function s = link2latex(g, phi, d, D)
% LINK2LATEX Generate LaTeX TikZ code from a directed graph with an endomap
%   s = link2latex(g, phi, d, D) returns a LaTeX string s.
%   - g is the complex coordinate vector (x + i*y) for each node
%   - phi is the endomap that defines the graph structure
%   - d is the domain map (node indices for the start of the edges)
%   - D is an optional vector of labels (defaults to node indices)

% Compute c from d and phi
c = d(phi);

% Handle optional labels
if nargin < 4
    D = 1:numel(g);
end

% Extract real and imaginary parts of coordinates
x = real(g);
y = imag(g);
n = numel(g);

% Header for LaTeX document
header = {
    '\begin{center}'
    '\begin{tikzpicture}[scale=1.5, every node/.style={circle,draw}, every path/.style={->, thick}]'
};

% Nodes
nodes = cell(n,1);
for i = 1:n
    nodes{i} = sprintf('\\node (n%d) at (%.2f, %.2f) {%d};', i, x(i), y(i), D(i));
end

% Edges (from d to c)
m = numel(d);
edges = cell(m,1);
for i = 1:m
    edges{i} = sprintf('\\path (n%d) edge (n%d);', d(i), c(i));
end

% Footer for LaTeX document
footer = {
    '\end{tikzpicture}'
    '\end{center}'
};

% Concatenate all lines into a single LaTeX string
lines = [header; nodes; edges; footer];
s = strjoin(lines, '\n');

% Optionally print to console
fprintf('%s\n', s);
end

%{
% Example
g = [1+1i, 2+2i, 3+1i, 2+3i]; % Example coordinates of nodes
phi = [2, 3, 4, 1]; % Endomap (mapping domain nodes to codomain nodes)
d = [1, 2, 3, 4]; % Domain map (node indices for the start of the edges)
D = 'ABCD'; % Optional labels for the nodes

latex_code = link2latex(g, phi, d, D);

%}
