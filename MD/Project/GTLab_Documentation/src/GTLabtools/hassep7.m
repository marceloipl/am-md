function [Sum] = hassep7(num)
% Topic: Number-Theorie
% Tag:   hasse diagramm
% Modification Log:
%   created:       25. Oct 2003 by A.Plattner (first version)
%   modified:      28. May 2010 by A.Plattner
% Despription:
% This function produces a hasse diagram
% Inputs:
% num: is the number you want to factorise
% Outputs:
% returns the sum of the factors, and the graph 
% Usage:
% - just run, the default num is 1008
% - hasse(num)
% Examples: ALIQUOT CYCLE, http://djm.cc/sociable.txt
% - hassep7(220)
% - K=14264, or K=14316 
% - K=hassep7(K)
% - K=14316, for h=1:28, K=hassep7(K), pause(1), end
% functions called:
% - 
% MATLAB: elaborated with R2007b, R2010a
format compact
clc, close all
if nargin == 0      % if run without parameter, take 1008 as num  
   num = 1008;
end
MEC=('MarkerEdgeColor');                    % define properties
MFC=('MarkerFaceColor');
MS =('MarkerSize');
W           = zeros(5,2);                   % make matrix W see example below
K2          = factor(num);
I           = 0;
while K2   ~= 0                             % produce table W with factors and exponend
	  I     = I+1;
	  A     = min(K2);
	  B     = max(find(K2 < (A+1)));
	  C     = find(K2 > A);
      K2    = K2(C);
      W(I,:)=[A B];
end
Xhoch=W(1,2)+1;     % dimension of X
Yhoch=W(2,2)+1;     % dimension of Y
Zhoch=W(3,2)+1;     % dimension of Z
X=(1:Xhoch);        % vector of X dimension
Y=(1:Yhoch);        % vector of Y dimension
Z=(1:Zhoch);        % vector of Z dimension
for m =Z                                    % produce the grid for the diagram
    Zx=ones(1,Xhoch);
    Zx=Zx*m;
    K=meshgrid(Y,X);
    for n=Y                                 % plot lines and speres in X-direction 
        plot3([X],K(:,n),Zx,'-r',[X],K(:,n),Zx,'o',...
              'LineWidth',1.0,MEC,'b',MFC,'b',MS,7); hold on
    end
    Zy=ones(1,Yhoch);
    Zy=Zy*m;
    K=meshgrid(X,Y);
    for n=X                                 % plot lines and speres in Y-direction
        plot3(K(:,n),[Y],Zy,'-r',K(:,n),[Y],Zy,'o',...
              'LineWidth',1.0,MEC,'b',MFC,'b',MS,7); hold on
    end
end
for m=Y
    Zz=ones(1,Zhoch);
    Zz=Zz*m;
    K=meshgrid(X,Z);
    for n=X                                 % plot lines in Z-direction
        plot3(K(:,n),Zz,[Z],'-r','LineWidth',1.0); hold on
    end
end
set(gcf,'Color',[1 1 1])                    % background color 'White'
view(-48,26);
axis equal
axis off
axis vis3d
xlim([0.5 W(1,2)+1.5])
Dim=length(find(W(:,1)));
r=0;                                        % labeling until dimension 3
b=0;                                        % factorise until dimension 5
for h=0:W(5,2)                              % dim = 5
    for k=0:W(4,2)                          % dim = 4
        for l=0:W(3,2)                      % dim z
            for m=0:W(2,2)                  % dim y
                for n=0:W(1,2)              % dim x
                    r = r+1;
                    a = ((W(5,1)^h)*(W(4,1)^k)*W(3,1)^l)*(W(2,1)^m)*(W(1,1)^n);
                    TR(r)=a;
                    if  length(find(W(:,1))) <= 3
                        text(n+1,m+.85,l+1,num2str(a),'FontWeight','bold','Color','k');
                    end
                end
            end
        end
    end
end    
% Print out to screen
fid=1;                                      % format output in command-window
fprintf(fid,'%12s','num = ')
fprintf(fid,'%2i\n\n',num)
fprintf(fid,'%12s','dimension = ')
fprintf(fid,'%2i\n\n',Dim)
% fprintf(fid,'%0i ',sort(TR))
  fprintf(fid,'%12s','factors = ')
fprintf(fid,'\n%12i %12i %12i %12i',sort(TR))
  display(char(13))
% W = W(1:length(find(W(:,1))),:)
TRS = num2str(sort(TR));
Sum = sum(TR)-num;
if Dim >= 4 
   disp('Warning: exceedes more then 3 dimesions')
   disp('Therefore, the Hasse-Graph is not complete and there is no labeling')
   disp('The factors and the Sum until dimension 5 are correct')
end
T(1,1)=strcat({'Hasse-Diagram of # '},{num2str(num)});
T(1,2)={TRS};
title(T,'FontSize',9,'Color','b')