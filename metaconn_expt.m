[nNodes,nTime] = size(data);


%% Calculate MTD (http://www.github.com/macshing/coupling/)
window = 10; % parameter that controls the window length of the simple moving average
[~,mtd] = coupling(data',window);


%% Multi-slice community detection (http://mucha.web.unc.edu/networks/)

A{nTime} = [];

for t = 1:nTime
  A{t} = mtd(:,:,t);
end

N = length(A{1});
T = length(A);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
S2 = zeros(nNodes,nTime);

% Temporal clustering analysis (requires code from: http://netwiki.amath.unc.edu/GenLouvain/GenLouvain)

gamma = 1; %parameter that defines the size of modules detected
omega = 1; %parameter that defines the strength of communities over time

for s=1:T
  k=sum(A{s});
  twom=sum(k);
  twomu=twomu+twom;
  indx=[1:N]+(s-1)*N;
  B(indx,indx)=A{s}-gamma*k'*k/twom;
end

twomu=twomu+2*omega*N*(T-1);
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B);
Q2 = Q/twomu;
S2(:,:) = reshape(S,N,T); %time-resolved community assignment


%% flexibility
%code from http://commdetect.weebly.com/uploads/4/9/5/9/49593677/flexibility.m

flex = flexibility(S2,'temp');




%% graph theoretical measures (requires code from: https://sites.google.com/site/bctnet/Home)

P = zeros(nNodes,nTime);
Z = zeros(nNodes,nTime);
q = zeros(nTime,1);

for t = 1:nTime
  P(:,t) = participation_coef_sign(sma(:,:,t),S2(:,t));
  Z(:,t) = module_degree_zscore(sma(:,:,t),S2(:,t));
end


%% cartographic profile

ybins = 8.5:-.105:-2;
xbins = 0:0.01:1.0;
xNumBins = numel(xbins); yNumBins = numel(ybins);
cart_profile = zeros(xNumBins,yNumBins,nTime);

for t = 1:nTime
  Xi = round( interp1(xbins, 1:xNumBins, P(:,t), 'linear', 'extrap') );
  Yi = round( interp1(ybins, 1:yNumBins, Z(:,t), 'linear', 'extrap') );
  Xi = max( min(Xi,xNumBins), 1);
  Yi = max( min(Yi,yNumBins), 1);
  cart_profile(:,:,t) = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);
end

