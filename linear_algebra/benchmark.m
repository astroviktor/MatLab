% Part b)

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
%     disp([' GE solution for system of size ',num2str(nlarge),' takes average time ',num2str(testtimesgauss(in)),' s']);
end %for

figure(1);
plot(nvals,testtimesgauss,'-bo','LineWidth',1.2)
xlabel('system size');
ylabel('time to solve (s)');
title('Empirically Determined Performance');

disp('Start of tests for Jacobi iteration');
tol=1e-9;
testtimesjacobi=zeros(size(nvals));
for in=1:numel(nvals)
    nlarge=nvals(in);
    Blarge=diag(-1*ones(nlarge-1,1),-1)+diag(-1*ones(nlarge-1,1),1)+diag(4*ones(nlarge,1),0);    %this must be diagonally dominant or else the method won't converge
    blarge=ones(nlarge,1);

    for irep=1:lrep     %benchmark will repeat the same solution several times to eliminate random variations from CPU load, etc.
        tstart=cputime;
        x0=randn(nlarge,1);
        [xit,iterations]=Jacobi(x0,Blarge,blarge,tol,false);
        tend=cputime;
        testtimesjacobi(in)=testtimesjacobi(in)+(tend-tstart)/lrep;
    end %for
%     disp([' JI solution for system of size ',num2str(nlarge),' takes average time ',num2str(testtimesjacobi(in)),' s']);
end %for

figure(1);
hold on
plot(nvals,testtimesjacobi,'-ro','LineWidth',1.2)
xlabel('system size');
ylabel('time to solve (s)');
legend('Gauss elim.','Jacobi it.')
title('Empirically Determined Performance');

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