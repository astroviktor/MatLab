  function [r, v, jd, coe, coe2] = PlanetData ...
    (planet_id, year, month, day, hour, minute, second, mu)

%   -- This code was adapted from planet_elements_and_sv.m --
% 
%   This function calculates the orbital elements and the state  
%   vector of a planet from the date (year, month, day)
%   and universal time (hour, minute, second).
%     
% 	INPUTS:
%   ==============================
%               
%   planet_id - planet identifier:
%                1 = Mercury
%                2 = Venus
%                3 = Earth
%                4 = Mars
%                5 = Jupiter
%                7 = Uranus
%                8 = Neptune
%                9 = Pluto
%   year      - range: 1901 - 2050
%   month     - range: 1 - 12
%   day       - range: 1 - 31
%   hour      - range: 0 - 23
%   minute    - range: 0 - 60
%   second    - range: 0 - 60
%   mu        - gravitational parameter of the sun (km^3/s^2)
%                          
%   OUTPUTS
%   ===============================
%   r         - heliocentric position vector
%   v         - heliocentric velocity vector
%   coe       - vector of MAIN heliocentric orbital elements
%               [a  incl  RA  w  e  TA  ] 
%               where
%                a     = semimajor axis                      (km)
%                incl  = inclination                         (rad)
%                RA    = right ascension                     (rad)
%                w     = argument of perihelion              (rad)
%                e     = eccentricity
%                TA    = true anomaly                        (rad)
%   coe2      - vector of ADDITIONAL heliocentric orbital elements
%                h     = angular momentum                    (km^2/s)
%                w_hat = longitude of perihelion ( = RA + w) (rad)
%                L     = mean longitude ( = w_hat + M)       (rad)
%                M     = mean anomaly                        (rad)
%                E     = eccentric anomaly                   (rad)
%    
%                 
%   Included subfunctions: J0, kepler_E, sv_from_coe, planetary_elements
%
% --------------------------------------------------------------------
%

%{         
  INTERNAL VARIABLES

  deg       - conversion factor between degrees and radians
  pi        - 3.1415926...

  j0        - Julian day number of the date at 0 hr UT
  ut        - universal time in fractions of a day
  jd        - julian day number of the date and time
 
  J2000_coe - row vector of J2000 orbital elements from Table 9.1
  rates     - row vector of Julian centennial rates from Table 9.1
  t0        - Julian centuries between J2000 and jd
  elements  - orbital elements at jd
%}

% Default value for mu of sun
if( nargin<8 )
  mu = 1.327e11; 
end

% By default - set hour, minute, second to 0 if not provided
if( nargin<7 )
  second = 0;
  if( nargin<6 )
    minute = 0;
    if( nargin<5 )
      hour = 0;
    end
  end
end

% FLAG - choose to force all orbits to be in ecliptic plane or not
ecliptic = 1;

% conversion constant
deg    = pi/180;

%...Equation 5.48:
j0     = J0(year, month, day);

ut     = (hour + minute/60 + second/3600)/24;

%...Equation 5.47
jd     = j0 + ut;

%...Obtain the data for the selected planet from Table 8.1:
[J2000_coe, rates] = planetary_elements(planet_id);

%...Equation 8.93a:
t0     = (jd - 2451545)/36525;

%...Equation 8.93b:
elements = J2000_coe + rates*t0;

a      = elements(1);
e      = elements(2);

%...Equation 2.71:
h      = sqrt(mu*a*(1 - e^2));

% Inclination - keep it? or force to be ecliptic?
if( ecliptic )
  incl = 0;
else
  incl = elements(3)*deg;
end

%...Reduce the angular elements to within the range 0 - 360 degrees:
RA     = mod(elements(4),360)*deg;
w_hat  = mod(elements(5),360)*deg;
L      = mod(elements(6),360)*deg;
w      = mod(w_hat - RA ,2*pi);
M      = mod(L - w_hat  ,2*pi);

%...Algorithm 3.1 (for which M must be in radians)
E      = kepler_E(e, M); %rad

%...Equation 3.13 (converting the result to degrees):
TA     = mod(2*atand(sqrt((1 + e)/(1 - e))*tan(E/2)),360)*deg;
  
coe    = [a incl RA w e TA];
coe2   = [h w_hat L M E];

%...Algorithm 4.5:
[r, v] = sv_from_coe([h e RA incl w TA], mu);
r = r(:);
v = v(:);
return

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function [J2000_coe, rates] = planetary_elements(planet_id)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
  This function extracts a planet's J2000 orbital elements and
  centennial rates from Table 8.1.
 
  planet_id      - 1 through 9, for Mercury through Pluto
 
  J2000_elements - 9 by 6 matrix of J2000 orbital elements for the nine
                   planets Mercury through Pluto. The columns of each 
                   row are:
                     a     = semimajor axis (AU)
                     e     = eccentricity
                     i     = inclination (degrees)
                     RA    = right ascension of the ascending
                             node (degrees)
                     w_hat = longitude of perihelion (degrees)
                     L     = mean longitude (degrees)
 
  cent_rates     - 9 by 6 matrix of the rates of change of the 
                   J2000_elements per Julian century (Cy). Using "dot"
                   for time derivative, the columns of each row are:
                     a_dot     (AU/Cy)
                     e_dot     (1/Cy)
                     i_dot     (deg/Cy)
                     RA_dot    (deg/Cy)
                     w_hat_dot (deg/Cy)
                     Ldot      (deg/Cy)
 
  J2000_coe      - row vector of J2000_elements corresponding
                   to "planet_id", with au converted to km
  rates          - row vector of cent_rates corresponding to
                   "planet_id", with au converted to km              
 
  au             - astronomical unit (149597871 km)
%}
% --------------------------------------------------------------------

%---- a --------- e -------- i -------- RA --------- w_hat ------- L ------

J2000_elements = ...
[0.38709927  0.20563593  7.00497902  48.33076593  77.45779628  252.25032350
 0.72333566  0.00677672  3.39467605  76.67984255 131.60246718  181.97909950 
 1.00000261  0.01671123 -0.00001531   0.0        102.93768193  100.46457166 
 1.52371034  0.09339410  1.84969142  49.55953891 -23.94362959 	-4.55343205
 5.20288700  0.04838624  1.30439695 100.47390909  14.72847983 	34.39644501
 9.53667594  0.05386179  2.48599187 113.66242448  92.59887831 	49.95424423
19.18916464  0.04725744  0.77263783  74.01692503 170.95427630  313.23810451
30.06992276  0.00859048  1.77004347 131.78422574  44.96476227  -55.12002969 
39.48211675  0.24882730 17.14001206 110.30393684 224.06891629  238.92903833];

cent_rates = ... 
[0.00000037  0.00001906 -0.00594749 -0.12534081  0.16047689  149472.67411175 
 0.00000390 -0.00004107 -0.00078890 -0.27769418  0.00268329	  58517.81538729  
 0.00000562 -0.00004392 -0.01294668  0.0         0.32327364   35999.37244981  
 0.0001847 	 0.00007882 -0.00813131 -0.29257343  0.44441088   19140.30268499  
-0.00011607 -0.00013253 -0.00183714  0.20469106	 0.21252668    3034.74612775 
-0.00125060 -0.00050991  0.00193609 -0.28867794 -0.41897216    1222.49362201
-0.00196176 -0.00004397 -0.00242939  0.04240589  0.40805281 	428.48202785 
 0.00026291  0.00005105  0.00035372 -0.00508664 -0.32241464 	218.45945325 
-0.00031596  0.00005170  0.00004818 -0.01183482 -0.04062942 	145.20780515]; 
 
J2000_coe      = J2000_elements(planet_id,:);
rates          = cent_rates(planet_id,:);

%...Convert from AU to km:
au             = 149597871; 
J2000_coe(1)   = J2000_coe(1)*au;
rates(1)       = rates(1)*au;

end %planetary_elements

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function E = kepler_E(e, M)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
  This function uses Newton's method to solve Kepler's 
  equation  E - e*sin(E) = M  for the eccentric anomaly,
  given the eccentricity and the mean anomaly.

  E  - eccentric anomaly (radians)
  e  - eccentricity, passed from the calling program
  M  - mean anomaly (radians), passed from the calling program
  pi - 3.1415926...

  User m-functions required: none
%}
% ----------------------------------------------

%...Set an error tolerance:
error = 1.e-8;

%...Select a starting value for E:
if M < pi
    E = M + e/2;
else
    E = M - e/2;
end

%...Iterate on Equation 3.17 until E is determined to within
%...the error tolerance:
ratio = 1;
while abs(ratio) > error
    ratio = (E - e*sin(E) - M)/(1 - e*cos(E));
    E = E - ratio;
end

end %kepler_E

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function [r, v] = sv_from_coe(coe,mu)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
  This function computes the state vector (r,v) from the
  classical orbital elements (coe).
 
  mu   - gravitational parameter (km^3;s^2)
  coe  - orbital elements [h e RA incl w TA]
         where
             h    = angular momentum (km^2/s)
             e    = eccentricity
             RA   = right ascension of the ascending node (rad)
             incl = inclination of the orbit (rad)
             w    = argument of perigee (rad)
             TA   = true anomaly (rad)
  R3_w - Rotation matrix about the z-axis through the angle w
  R1_i - Rotation matrix about the x-axis through the angle i
  R3_W - Rotation matrix about the z-axis through the angle RA
  Q_pX - Matrix of the transformation from perifocal to geocentric 
         equatorial frame
  rp   - position vector in the perifocal frame (km)
  vp   - velocity vector in the perifocal frame (km/s)
  r    - position vector in the geocentric equatorial frame (km)
  v    - velocity vector in the geocentric equatorial frame (km/s)

  User M-functions required: none
%}
% ----------------------------------------------

h    = coe(1);
e    = coe(2);
RA   = coe(3);
incl = coe(4);
w    = coe(5);
TA   = coe(6);

%...Equations 4.45 and 4.46 (rp and vp are column vectors):
rp = (h^2/mu) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);
vp = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);

%...Equation 4.34:
R3_W = [ cos(RA)  sin(RA)  0
        -sin(RA)  cos(RA)  0
            0        0     1];

%...Equation 4.32:
R1_i = [1       0          0
        0   cos(incl)  sin(incl)
        0  -sin(incl)  cos(incl)];

%...Equation 4.34:
R3_w = [ cos(w)  sin(w)  0 
        -sin(w)  cos(w)  0
           0       0     1];

%...Equation 4.49:
Q_pX = (R3_w*R1_i*R3_W)';

%...Equations 4.51 (r and v are column vectors):
r = Q_pX*rp;
v = Q_pX*vp;

%...Convert r and v into row vectors:
r = r';
v = v';

end
% sv_from_coe

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function j0 = J0(year, month, day)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%{
  This function computes the Julian day number at 0 UT for any year
  between 1900 and 2100 using Equation 5.48.
 
  j0    - Julian day at 0 hr UT (Universal Time)
  year  - range: 1901 - 2099
  month - range: 1 - 12
  day   - range: 1 - 31
 
  User m-functions required: none
%}
% ----------------------------------

j0 = 367*year - fix(7*(year + fix((month + 9)/12))/4) ...
    + fix(275*month/9) + day + 1721013.5;
   		
end %J0

end %planet_elements_and_sv
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
