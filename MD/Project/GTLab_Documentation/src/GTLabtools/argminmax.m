function y=argminmax(x,str)
% y=argminmax(x,str)
% computed [~,i]=sort(x) and returns:
% (a) i(1) if str=='min'; 
% (b) i(end) if str=='max'

[~,i]=sort(x);
switch str
    case 'min'
        y=i(1,:);
    case 'max'
        y=i(end,:);
end
end % argminmax