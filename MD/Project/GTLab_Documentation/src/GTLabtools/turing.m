function Out=turing(tape,head,Write,Move,Transit,varargin)
% Out=turing(tape,head,Write,Move,Transit)
% Fonte: Matematica Discreta 23 - PL (Teams)
%
% simulates the behaviour of a Turing Machine with initial input as in
% 'tape', the input 'head' indicates the initial position and the three 
% matrices of the same size Write, Move and Transit define the quintuples 
% q=(si,read,write,move,sj)
% with the usual meaning that when at the internal state si, if reading
% 'read' then writes 'write', moves by 'move' and transits to the new state
% sj
% the three matrices Write, Move, Transit have to be of the same size, and
% the lines are indexed with the alphabet letters [0, 1, 2, ...] while the
% states are indexes by the columns as [1,2,3,...] with NaN as the halting
% state, the initial state is always assumed to be 1
%
% Example:
%{
tape=zeros(1,10); head=4; Write=[1 0 1; 1 1 1]; Move=[1 1 -1; 1 1 -1];
Transit=[2 3 3; NaN 2 1]; % The busy beaver TM with three states
%}
% Use the syntax
% Out=turing(zeros(1,10),4,Write,Move,Transit,'disp')
% for a display of the output at each step, needs the function mydigraph.m
% and it is limited to 16 states

version='27May2025'
%version='20May2025'
%version 'Matematica Discreta 23 - PL (Teams)'

switch version
case '27May2025'
%{
%example 
tape=20; head=10;
Write=[1 0 1; 1 1 1]; Move=[1 1 -1; 1 1 -1];
Transit=[2 3 3; 0 2 1]; % The busy beaver TM with three states 0=halt
Out=turing(20,10,Write,Move,Transit)
%}
if numel(tape)==1
    tape=sparse(1,ceil(abs(tape)));
end
Phi=(Write+1).*Move + 1i*Transit;
z=0+1i % initial state is 1
stop=0;
while ~stop
    z=Phi(tape(head)+1,imag(z));
    tape(head)=abs(real(z))-1;
    head=head+sign(real(z));
    if imag(z)==0||(tape(head)<0)
        stop=1;
    end    
end
Out=tape;


case '20May2025'

%{
%example 
tape=zeros(1,10); head=4; Write=[1 0 1; 1 1 1]; Move=[1 1 -1; 1 1 -1];
Transit=[2 3 3; 0 2 1]; % The busy beaver TM with three states
%}

Phi=(Write+1).*Move + 1i*Transit;
tape=sparse(1,100);
state=1;
head=50;
stop=0;
while ~stop
    z=Phi(tape(head)+1,state);
    tape(head)=abs(real(z))-1;
    state=imag(z);
    stop=isequal(state,0)||(tape(head)<0);
    head=head+sign(real(z));
end
Out=tape;


case 'Matematica Discreta 23 - PL (Teams)'
if nargin==5
Out=tape;
state=1;
while ~isnan(state)
    read=Out(head);
    Out(head)=Write(read+1,state);
    head=head+Move(read+1,state);
    state=Transit(read+1,state);
end

else % in the case that the sixth input argument is varargin=='disp'
nlc=size(Write);
if ~isequal(size(Move),nlc) || ~isequal(size(Transit),nlc)
    error('The input matrices Write, Move and Transit must have the same size')
end
%[s,t,u]=mydigraph(Write(:),tape(:));
[s,t,u]=digraph(Write(:),tape(:)); %19May2025
W=reshape(s,size(Write));
T=t';

state=1;

while ~isnan(state)

    N=nan(size(T)); N(head)=dec2hex(state); %limited to 16 states
    disp(char([T+64;N]))

    read=T(head);
    T(head)=W(read,state);
    head=head+Move(read,state);
    state=Transit(read,state);
end
disp(char(T+64))
Out=u(T)';
end

end % Switch case
