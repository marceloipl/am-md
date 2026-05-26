function gp=slerp(a,b,t)
  % gp=slerp(a,b,t)
  %
  % given two points a and b on a sphere and a parameter t in [0,1], returns the
 %  geometrical path gp connecting a and b at instance t along the sphere's surface

 % compute the unitary dot product of a and b
 udotab=dot(a,b)./(norm(a).*norm(b));

 % adjusting the inputs
 a=a(:)'; b=b(:)'; t=t(:);

 % analyse special cases first
 if udotab==1  % a and b are parallel
   pg=a+t*(b-a);
 elseif udotab==-1 % antipodal
   disp('the case of antipodals is not implemented yet')
   gp=[];
 else
   omega= acos(udotab);
   gp= sin((1-t)*omega)/sin(omega)*a + sin(t*omega)/sin(omega)*b;
 end
