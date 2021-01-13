      clear

      i4digit = input('What is the 4 digit number?')
      ipts    = input('Enter the number of points on each side of airfoil')

      icmax =floor(i4digit/1000.0);
      ixcmax=floor(i4digit/100.0)-icmax*10;
      itmax =i4digit-ixcmax*100-icmax*1000;

      cmax = icmax*0.01;
      xcmax=ixcmax*0.1;
      tmax = itmax*0.01;

      for i=1:ipts;
      theta = pi*(i-1)/(ipts-1);
      x(i) = 0.5*(1.0-cos(theta));
      yt(i) = tmax/0.20*(0.29690*sqrt(x(i))-0.12600*x(i)...
             -0.35160*x(i)*x(i)+0.28430*x(i)^3-0.10360*x(i)^4);
      yc(i) = 0.0;
      if (icmax~=0);
       if (x(i)<=xcmax);
        yc(i) = cmax/(xcmax*xcmax)*(2.0*xcmax*x(i)-x(i)*x(i));
       else;
        yc(i) = cmax/(1.0-xcmax)^2*((1.0-2.0*xcmax)+2.0*xcmax*x(i)-x(i)*x(i));
       end;
      end;
      end;

      mp1 = ipts*2-1;
      ii = 1;
      for i=ipts:-1:1;
      xb(ii) = x(i);
      yb(ii) = yc(i)-yt(i);
      ii = ii+1;
      end;
      for i=2:ipts;
      xb(ii) = x(i);
      yb(ii) = yc(i)+yt(i);
      ii = ii+1;
      end;

      plot(xb,yb);

      xb = xb';
      yb = yb';

      m = mp1-1;

      alpha = 8.0*pi/180;

      for i = 1: m;
      ip1  = i + 1;
      x(i) = 0.5*(xb(i)+xb(ip1));
      y(i) = 0.5*(yb(i)+yb(ip1));
      s(i) = sqrt( (xb(ip1)-xb(i))^2 + (yb(ip1)-yb(i))^2 );
      theta(i)  = atan2( (yb(ip1)-yb(i)), (xb(ip1)-xb(i)) );
      sine(i)   = sin( theta(i) );
      cosine(i) = cos( theta(i) );
      rhs(i)    = sin( theta(i)-alpha );
      end;

      for i = 1: m;
      for j = 1: m;
      if ( i ~= j );
      a = -(x(i)-xb(j))*cosine(j) - (y(i)-yb(j))*sine(j);
      b = (x(i)-xb(j))^2 + (y(i)-yb(j))^2;
      c = sin( theta(i)-theta(j) );
      d = cos( theta(i)-theta(j) );
      e = (x(i)-xb(j))*sine(j) - (y(i)-yb(j))*cosine(j);
      f = log( 1.0 + s(j)*(s(j)+2.*a)/b );
      g = atan2( e*s(j), b+a*s(j) );
      p = (x(i)-xb(j)) * sin( theta(i)-2.*theta(j) )...
        + (y(i)-yb(j)) * cos( theta(i)-2.*theta(j) );
      q = (x(i)-xb(j)) * cos( theta(i)-2.*theta(j) )...
        - (y(i)-yb(j)) * sin( theta(i)-2.*theta(j) );
      cn2(i,j) = d + .5*q*f/s(j) - (a*c+d*e)*g/s(j);
      cn1(i,j) = .5*d*f + c*g - cn2(i,j);
      ct2(i,j) = c + .5*p*f/s(j) + (a*d-c*e)*g/s(j);
      ct1(i,j) = .5*c*f - d*g - ct2(i,j);
      else;
      cn1(i,j) = -1.0;
      cn2(i,j) =  1.0;
      ct1(i,j) =  0.5*pi;
      ct2(i,j) =  0.5*pi;
      end;
      end;
      end;

      for i = 1: m;
      an(i,1)   = cn1(i,1);
      an(i,mp1) = cn2(i,m);
      at(i,1)   = ct1(i,1);
      at(i,mp1) = ct2(i,m);
      for j = 2: m;
      an(i,j) = cn1(i,j) + cn2(i,j-1);
      at(i,j) = ct1(i,j) + ct2(i,j-1);
      end;
      end;

      an(mp1,1)   = 1.0;
      an(mp1,mp1) = 1.0;

      for j = 2: m;
      an(mp1,j) = 0.0;
      end;
      rhs(mp1)  = 0.0;

      rhs = rhs';
      gama = an\rhs;

      for i = 1: m;
      v(i) = cos( theta(i)-alpha );
      for j = 1: mp1;
      v(i) = v(i) + at(i,j)*gama(j);
      cp(i) = 1.0 - v(i)^2;
      end;
      end;

      plot(x,-cp)
