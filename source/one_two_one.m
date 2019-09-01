function sv = one_two_one(v)

% apply a 1-2-1 filter in one dimension

  %error(nargchk(1, 1, nargin));

  %v         = v(:);
  sv        = zeros(size(v));
  nv        = length(v);

  if nv<3
   disp('too few points for filter implementation, returning the same thing')
   sv = v;
   return
  end

  for j=2:nv-1
   sv(j) = 0.5*v(j)+0.25*v(j+1)+0.25*v(j-1);
  end

  sv(1) = 0.75*v(1)+0.25*v(2);
  sv(nv) = 0.75*v(nv)+0.25*v(nv-1);

