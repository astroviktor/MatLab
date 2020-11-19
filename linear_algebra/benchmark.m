<<<<<<< HEAD
% Part b)
=======
% Evaluate performance and scaling of Gaussian elimination and Jacobi iteration
%    by solving systems of different size and timing the solves
>>>>>>> 2782018b8d64a2ee8ff310ac767e4e9c7df7c489

% Evaluate performance and scaling of Gaussian elimination, Jacobi
% iteration and Thomas algorithm
%    by solving systems of different size and timing the solves
clear all;close all;clc
nvals=1:1:150;
testtimesgauss=zeros(size(nvals));
lrep=10;     %how many times to repeat each test

disp('Start of tests of Gaussian-elimination scaling');
for in=1:numel(nvals)
    nlarge=nvals(in);
    Blarge=diag(-1*ones(nlarge-1,1),-1)+diag(-1*ones(nlarge-1,1),1)+diag(4*ones(nlarge,1),0);    %this must be diagonally dominant or else the method won't converge
    blarge=ones(nlarge,1);
    
    for irep=1:lrep     %benchmark will repeat the same solution several times to eliminate random variations from CPU load, etc.
        tstart=cputime;
        [Blargemod,ordlarge]=Gauss_elim(Blarge,blarge);
        xlarge=backsub(Blargemod(ordlarge,:));
        tend=cputime;
        testtimesgauss(in)=testtimesgauss(in)+(tend-tstart)/lrep;
    end %for
<<<<<<< HEAD
%     disp([' GE solution for system of size ',num2str(nlarge),' takes average time ',num2str(testtimesgauss(in)),' s']);
end %for

figure(1);
plot(nvals,testtimesgauss,'-bo','LineWidth',1.2)
=======
    disp([' GE solution for system of size ',num2str(nlarge),' takes average time ',num2str(testtimes(in)),' s']);
end %for

figure(1);
plot(nvals,testtimes,'o','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','blue')
>>>>>>> 2782018b8d64a2ee8ff310ac767e4e9c7df7c489
xlabel('system size');
ylabel('time to solve (s)');
title('Empirically Determined Performance');

disp('Start of tests for Jacobi iteration');
tol=1e-9;
<<<<<<< HEAD
testtimesjacobi=zeros(size(nvals));
=======
testtimes=zeros(size(nvals));
>>>>>>> 2782018b8d64a2ee8ff310ac767e4e9c7df7c489
for in=1:numel(nvals)
    nlarge=nvals(in);
    Blarge=diag(-1*ones(nlarge-1,1),-1)+diag(-1*ones(nlarge-1,1),1)+diag(4*ones(nlarge,1),0);    %this must be diagonally dominant or else the method won't converge
    blarge=ones(nlarge,1);

    for irep=1:lrep     %benchmark will repeat the same solution several times to eliminate random variations from CPU load, etc.
        tstart=cputime;
        x0=randn(nlarge,1);
        [xit,iterations]=Jacobi(x0,Blarge,blarge,tol,false);
        tend=cputime;
<<<<<<< HEAD
        testtimesjacobi(in)=testtimesjacobi(in)+(tend-tstart)/lrep;
    end %for
%     disp([' JI solution for system of size ',num2str(nlarge),' takes average time ',num2str(testtimesjacobi(in)),' s']);
=======
        testtimes(in)=testtimes(in)+(tend-tstart)/lrep;
    end %for
    disp([' JI solution for system of size ',num2str(nlarge),' takes average time ',num2str(testtimes(in)),' s']);
>>>>>>> 2782018b8d64a2ee8ff310ac767e4e9c7df7c489
end %for

figure(1);
hold on
<<<<<<< HEAD
plot(nvals,testtimesjacobi,'-ro','LineWidth',1.2)
=======
plot(nvals,testtimes,'^','LineWidth',2,'MarkerSize',20,'MarkerFaceColor','blue')
>>>>>>> 2782018b8d64a2ee8ff310ac767e4e9c7df7c489
xlabel('system size');
ylabel('time to solve (s)');
legend('Gauss elim.','Jacobi it.')
title('Empirically Determined Performance');
<<<<<<< HEAD

testtimesthomas=zeros(size(nvals));
disp('Start of Thomas method (tridiagonal)');
for in=1:numel(nvals)
    nlarge=nvals(in);
    Blarge=diag(-1*ones(nlarge-1,1),-1)+diag(-1*ones(nlarge-1,1),1)+diag(4*ones(nlarge,1),0);    %this must be diagonally dominant or else the method won't converge
    blarge=ones(nlarge,1);
    
    for irep=1:lrep     %benchmark will repeat the same solution several times to eliminate random variations from CPU load, etc.
        tstart=cputime;
        x=tridiag(Blarge,blarge);
        tend=cputime;
        testtimesthomas(in)=testtimesthomas(in)+(tend-tstart)/lrep;
    end %for
%     disp([' Thomas solution for system of size ',num2str(nlarge),' takes average time ',num2str(testtimesthomas(in)),' s']);
end %for

for i=1:length(testtimesgauss)
    if testtimesgauss(i)<testtimesjacobi(i)
       fprintf('For systems of size %i the Gaussian elimination is faster than the Jacobi iteration\n',i)
    end
end

figure(1);
hold on
plot(nvals,testtimesthomas,'-ko','LineWidth',1.2)
xlabel('system size');
ylabel('time to solve (s)');
legend('Gauss elim.','Jacobi it.','Thomas alg.')
title('Empirically Determined Performance');
=======
>>>>>>> 2782018b8d64a2ee8ff310ac767e4e9c7df7c489
