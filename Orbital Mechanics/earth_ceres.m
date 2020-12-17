%% Earth's Relative Motion about Ceres

close all
clear all

 image_file = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
  set(gca, 'NextPlot','add', 'Visible','off');
  
  axis equal;
  axis auto;
  axis vis3d;
  hold on;
  
  Re = 6378.14; % equatorial radius (km)
  [x, y, z] = ellipsoid(0, 0, 0, Re, Re, Re, 180);
  globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
  
  cdata = imread(image_file);
  set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');

   image_file2 = 'https://en.wikipedia.org/wiki/Ceres_(dwarf_planet)#/media/File:Ceres_-_RC3_-_Haulani_Crater_(22381131691)_(cropped).jpg';
  set(gca, 'NextPlot','add', 'Visible','off');
  
  axis equal;
  axis auto;
  axis vis3d;
  hold on;
  
  Rc = 469.73; % equatorial radius (km)
  [x, y, z] = ellipsoid(0, 0, 0, Rc, Rc, Rc, 180);
  globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
  
  cdata = imread(image_file2);
  set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');