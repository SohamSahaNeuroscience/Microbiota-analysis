[coeff, score, latent, explained, tsquared] = princomp(allmean);
figure;
scatter(coeff(:,1), coeff(:,2), 'filled', 'g'); hold on
scatter(coeff(:,2), coeff(:,3),'filled',   'r'); hold on
scatter(coeff(:,3), coeff(:,4),'filled', 'b'); hold on
hold off
a = rotatefactors(coeff(:,1:2),'Method','equamax');
figure;
biplot(a);


clearvars -except Index

% donor data mean
%d1 = mean(Index(1).FAMILY(:,1:4),2);%control T0
d2 = mean(Index(1).GENUS(:,9:13),2);%donor control T7
%d3 = mean(Index(1).FAMILY(:,5:8),2);%donor UCMS T0
d4 = mean(Index(1).GENUS(:,14:18),2);%donor UCMS T7
donor_mean = [d2 d4];

%receptor data mean
%r1 = mean(Index(1).FAMILY(:,19:20),2);%receptor T0
r2 = mean(Index(1).GENUS(:,32:36),2);%receptor T7
%r3 = mean(Index(1).FAMILY(:,21:23),2);%receptor T0 UCMS
r4 = mean(Index(1).GENUS(:,37:41),2);%receptor UCMS T7
receptor_mean = [r2 r4];

allmean = [donor_mean  receptor_mean];
vari = var(allmean');
[s in] = sort(vari, 'descend');
donor_mean_sorted = donor_mean(in,:);
receptor_mean_sorted = receptor_mean(in,:);

load('/Users/sohamsaha/Documents/QIIME/Eleni/RWBmap2.mat');
allmean_sorted = [donor_mean_sorted  receptor_mean_sorted];
imagesc(allmean_sorted);colormap(RWBmap);
caxis([0 0.25]);
[coeff, score] = pca(allmean_sorted');
figure;
scatter(coeff(:,1), coeff(:,2), 'filled', 'g','Linewidth', 12); 


%% without mean

%dd1 = Index(1).FAMILY(:,1:4);%control T0
dd2 = Index(1).GENUS(:,9:13);%donor control T7
%dd3 = Index(1).FAMILY(:,5:8);%donor UCMS T0
dd4 = Index(1).GENUS(:,14:18);%donor UCMS T7
donor_all = [dd2 dd4];

%receptor data mean
%rr1 = Index(1).FAMILY(:,19:20);%receptor T0
rr2 = Index(1).GENUS(:,32:36);%receptor T7
%rr3 = Index(1).FAMILY(:,21:23);%receptor T0 UCMS
rr4 = Index(1).GENUS(:,37:41);%receptor UCMS T7
receptor_all = [rr2 rr4];

alldats = [donor_all  receptor_all];
varia = var(alldats');
[ss ind] = sort(varia, 'descend');
donor_all_sorted = donor_all(ind,:);
receptor_all_sorted = receptor_all(ind,:);
alldats_sorted = [donor_all_sorted  receptor_all_sorted];
imagesc(alldats_sorted);colormap(RWBmap);
caxis([0 0.7]);

[coeff1, score1] = pca(alldats_sorted');
figure;
scatter(coeff1(1:5,1), coeff1(1:5,2), 'filled', 'g','Linewidth', 2, 'SizeData', 100); xlim([-1.6 1.6]);ylim([-1.6 1.6]); hold on
scatter(coeff1(6:10,1), coeff1(6:10,2), 'filled', 'r','Linewidth', 2, 'SizeData', 100);xlim([-1.6 1.6]);ylim([-1.6 1.6]);hold on
scatter(coeff1(11:15,1), coeff1(11:15,2), 'filled', 'b','Linewidth', 2, 'SizeData', 100);xlim([-1.6 1.6]);ylim([-1.6 1.6]); hold on
scatter(coeff1(16:20,1), coeff1(16:20,2), 'filled', 'k','Linewidth', 2, 'SizeData', 100);xlim([-1.6 1.6]);ylim([-1.6 1.6]);
hold off

% figure;
% scatter(coeff1(19:21,1), coeff1(19:21,2), 'filled', 'g','Linewidth', 2, 'SizeData', 100); xlim([-1.6 1.6]);ylim([-1.6 1.6]);hold on
% scatter(coeff1(22:27,1), coeff1(22:27,2), 'filled', 'r','Linewidth', 2, 'SizeData', 100); xlim([-1.6 1.6]);ylim([-1.6 1.6]);hold on
% scatter(coeff1(28:30,1), coeff1(28:30,2), 'filled', 'b','Linewidth', 2, 'SizeData', 100); xlim([-1.6 1.6]);ylim([-1.6 1.6]);hold on
% scatter(coeff1(31:34,1), coeff1(31:34,2), 'filled', 'k','Linewidth', 2, 'SizeData', 100);  xlim([-1.6 1.6]);ylim([-1.6 1.6]);
% hold off

%% euclidean dist




data = allmean_sorted;
ln_data = real(log(data));
blah = -(data .* ln_data);

for i = 1:size(blah,1)
    for j = 1:size(blah,2)
        if isnan(blah(i,j))== 1
            blah(i,j) = 0;
        end
    end
end
alpha = 1./(blah);
% 

figure;imagesc(donor_mean_sorted);colormap(RWBmap);caxis([-0.1 0.5]);
figure;imagesc(receptor_mean_sorted);colormap(RWBmap);caxis([-0.1 0.5]);


%%

v = var(x');
[ind, m] = sort(v,'descend');
x1 = x(m,:);

y = [x1(:,1:5) x1(:,11:15)];%CT
q = sum(y,2);
ind1 = find(q == 0);

y1 = [x1(:,6:10) x1(:,16:end)];%UCMS
q = sum(y1,2);
ind2 = find(q == 0);
IN = union(ind1, ind2);

y(IN,:) = [];
y1(IN,:) = [];

s = corr(y');
s1 = corr(y1');

figure;heatmap(s);
caxis([-0.5 1])
figure;heatmap(s1);
caxis([-0.5 1])