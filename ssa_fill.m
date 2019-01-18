function [x,R]=ssa_fill(x,M)
% [R]=ssa_fill(x,M)
% singular spectrum analysis fill in missing data algorithm 
% algorithm: schoellhammer 2001, geophysical research letters
% sanjay manohar 2016
ALLOWED_NANS = 0.7;                 % fraction of window that can be missing
mx=nanmean(x); x=x-mx;              % subtract mean
if isrow(x), x=x'; end              % ensure data is a column
N=size(x,1);                        % num data points
if ~exist('M','var'), M=floor(N/100); end % window size
%%%% begin SSA 
for i=1:M                           % for each delay, autocorrelate
  c(i)=nanmean(x(1:N-i+1).*x(i:N)); % simply ignore nans
end
[E,~] = eig( toeplitz(c) );
A=zeros(N,M);
for i=1:N-M+1;                      % project data onto eigs to get pc
  tmp    = x(i:i+M-1);              % get M-sized window to project
  if mean(isnan(tmp)) > ALLOWED_NANS
    A(i,:) = nan;                   % if there aren't enough data, pc=nan
  else
    ratio  = M/sum(~isnan(tmp));    % if nans, use proportionally
    e=E; % e(:,isnan(tmp))=0;
    tmp(isnan(tmp)) = 0;            % missing data contributes 0 to pc
    A(i,:) = tmp'*e *ratio;         % but scale up pc to make up for missing data
  end
end
A(N-M+2:N,:)=0; % pad bottom end with M-1 zeros
for k=1:M       % reconstruct from pc and eigs
  R(:,k)=filter(E(:,k),M,A(:,k));
end
for i=1:M-1     % Adjust first M-1 rows and last M-1 rows
  R(i,:)=R(i,:)*(M/i);
  R(N-i+1,:)=R(N-i+1,:)*(M/i);
end
R=fliplr(R);    % put most significant eigenvalue in front
x = sum(R,2);   % replace with fitted value
x = x + mx;     % add mean back
