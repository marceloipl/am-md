function [g,phi]=linkfigs(Fig,w)
% 14Oct2024 changes output to [phi,g]
  % [phi,g]=linkfigures(Fig,w)
  %
  % Returns a link structure for figures in the database, currently letter A and Number 1
  % added a star '*'

  switch Fig
    case 'A'
      phi=[2:8,1];
      g=[-1; 0; w; w-1; w/sqrt(2)-1; -conj(w)/sqrt(2); -conj(w); -conj(w)-1];
    case '1'
    phi=[2:8,1]; b=imag(w);
    g=[-1; 0; w; b*1i; pi*b*1i; pi*b*1i + w/pi; -1+pi*b*1i-conj(w)/pi; -1+pi*b*1i ];
  case '*'
      n=ceil(real(w)); m=ceil(imag(w));
      phi=reshape((1:n*m),[n m]);
      phi(1,:)=[]; phi(end+1,:)=phi(end,:); phi=phi(:);
      [a,b]=ndgrid(log((1:n)/n),(1:m)/m*2*pi);
      t=complex(a,b);
      g=exp(t(:));
    case 'circles'
      n=ceil(real(w)); m=ceil(imag(w))+1;
      phi=reshape((1:n*m),[n m]);
      phi=reshape(phi(:,[2:end 1]),[m*n 1]);
      [a,b]=ndgrid(log((1:n)/n),(1:m)/m*2*pi);
      t=complex(a,b);
      g=exp(t(:));
      case 'M'
      phi=[2 3 4 5 5];
      g=[0+0i, 0+1i, 0.5+0.5i, 1+1i, 1+0i];
    case 'V'
      phi=[2 3 4 5 5];
      g=[0+1i, 0+0.5i, 0.5+0i, 1+0.5i, 1+1i];
    case 'line'  % w=5+5i
      n=abs(ceil(real(w))); m=abs(ceil(imag(w)));
      if m>n, m=n; end
      if m<1, m=1; end
      phi=[2:n, m];
      g=linspace(0,1,n);
    otherwise
      disp([Fig ' is not a valid input string'])
      disp('The current possibilities are: A, 1, *, circles')
      phi=[]; g=[];
  end
