rng('default');  % For reproducibility
clusterTest_x = [gallery('uniformdata',[10 3],12); ...
    gallery('uniformdata',[10 3],13)+1.2; ...
    gallery('uniformdata',[10 3],14)+2.5];
clusterTest_y = [ones(10,1);2*(ones(10,1));3*(ones(10,1))]; % Actual classes

scatter3(clusterTest_x(:,1),clusterTest_x(:,2),clusterTest_x(:,3),100,clusterTest_y,'filled')
title('Randomly Generated Data in Three Clusters');

T = clusterdata(clusterTest_x, 'Linkage','ward','SaveMemory','on','Maxclust',4);

scatter3(clusterTest_x(:,1),clusterTest_x(:,2),clusterTest_x(:,3),100,T,'filled')

%%
X = rand(2000,3);

% rows = observations; columns = categories/dimensions
T = clusterdata(X,'Linkage','average','SaveMemory','off','Maxclust',4);

% scatter3(X(:,1),X(:,2),X(:,3),10,T,'filled')


%%
y = pdist(X);
% y_sq = squareform(y);

z = linkage(y);

clf, figure(1)
dendrogram(z)

%%
% z_pd = pdist(X, 'cityblock');
% z = linkage(z_pd, 'average');
z_pd = pdist(X, 'euclidean');
z = linkage(z_pd, 'ward');

figure(2)
dendrogram(z,'Reorder',order)

%%
z_pd = pdist(X, 'cityblock');
tree = linkage(z_pd, 'average');

figure(3)
[~,T_d] = dendrogram(tree, 4);

