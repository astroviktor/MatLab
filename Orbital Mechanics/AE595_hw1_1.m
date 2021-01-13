%ORBITAL MECHANICS
%Homework 2
%Problem 2 and 3

%picking random position and velocity vectors
r=[3482;1200;4392]; %km
v=[4;-2;6]; %km/s

%gravitational parameter for Earth
mu=398600; %Km^3/s^2

[a,i,omega,w,e,th,h] = OrbitalElementsFromRV( r, v, mu );

fprintf('The keplerian orbital elements for such position and velocity are:\nSemi mayor axis=%1.4d km\nInclination=%1.2d rad\nRAAN=%1.2d rad\nArgument of perigee=%.1i rad\nEccentricity=%i\nTrue anomaly=%1.2d rad\n\n',a,i,omega,w,e,th);


%Let's double check our work with getting r and v from the orbital elements

[r,v] = RVFromCOE( a,i,omega,w,e,th, mu );

fmt = ['The velocity vector is confirmed to be: [', repmat('%g, ', 1, numel(v)-1), '%g] km/s\n'];
fprintf(fmt, v)
fmt = ['The position vector is confirmed to be: [', repmat('%g, ', 1, numel(r)-1), '%g] km\n'];
fprintf(fmt, r)


