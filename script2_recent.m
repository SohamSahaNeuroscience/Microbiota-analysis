close all
clear all
clc

rootdir = '/Users/sophiesoham/Desktop/Eleni/ALL_ANALYSIS/Matlab variables/Family';
cd(rootdir)

load('data_final.mat');
a =log2(datapoints);
[x y] = size(a);
for i = 1:x
    for j = 1:y
        if a(i,j) == -Inf
            a(i,j) = 0;
        end
    end
end

sum_data = sum(a,2);
index = find(sum_data == 0);
a(index,:)=[];
corrected_data = a(:,1:42);
clear a
donor = corrected_data(:,1:18);
% donor data mean
d1 = donor(:,1:4);
d2 = donor(:,5:8);
d3 = donor(:,9:13);
d4 = donor(:,14:18);
donor_mean = [d1 d2 d3 d4];


receptor = corrected_data(:,19:end);
%receptor data mean
r1 = receptor(:,7:10);
r2 = receptor(:,11:14);
r3 = receptor(:,15:19);
r4 = receptor(:,20:24);
receptor_mean = [r1 r2 r3 r4];

save('corrected_data.mat','corrected_data');
save('donor.mat','donor');
save('receptor.mat','receptor');
save('donor_mean.mat','donor_mean');
save('receptor_mean.mat','receptor_mean');


clear i j x y


x = d3;
y = d4;
x1 = r3;
x2 = r4;
norm_donor = [x y];
norm_receptor = [x1 x2];
total = [norm_donor norm_receptor];
for i = 1:size(total,2)
    for j = 1:size(total,1)
        if total(j,i)<0
            total(j,i)=-1;
        end
    end
end
clear d1 d2 d3 d4
clear r1 r2 r3 r4

total_round = round(total);
clear x y x1 x2 i j
save('difference.mat','total_round');

%% for indexing bacterial OTUids

load('family.mat');
x = eval('Family');
x(index) = [];%index = from the prev indexing of data

load('otuid.mat');
y = eval('OTUID');
y(index) = [];
bacteria = {'Actinobacteria','Bacteroidetes','Chlamydiae','Cyanobacteria','Deferribacteres','Firmicutes','Planctomycetes','Proteobacteria','Saccharibacteria','Tenericutes','Verrucomicrobia'};
bacteria = eval('bacteria');

for i = 1:length(bacteria)
Index(i).names = find(strcmp(y, bacteria(i)));
Index(i).fam = x(Index(i).names);
end
save('Index.mat','Index');
clear i OTUID

for i = 1:length(Index)
    Index(i).donor_mean = donor_mean(Index(i).names,:);
    Index(i).receptor_mean = receptor_mean(Index(i).names,:);
end
clear i
for i = 1:length(Index)
    Index(i).donor = donor_mean(Index(i).names,:);
    Index(i).receptor = receptor_mean(Index(i).names,:);
end

blah = [];
blah1 = [];
for i= 1:length(Index)
    blah = [blah; Index(i).donor_mean];
    Donor_indexed = blah;
    blah1 = [blah1; Index(i).receptor_mean];
    Receptor_indexed = blah1;
end
clear blah blah1

blah = [];
blah1 = [];
for i= 1:length(Index)
    blah = [blah; Index(i).donor];
    Donor = blah;
    blah1 = [blah1; Index(i).receptor];
    Receptor = blah1;
end
clear blah blah1

save('Donor_indexed.mat','Donor_indexed');
save('Receptor_indexed.mat','Receptor_indexed');
save('Index.mat','Index');

x1 = Donor_indexed(:,3)-Donor_indexed(:,1);
y1 = Donor_indexed(:,4)-Donor_indexed(:,2);
x2 = Receptor_indexed(:,3)-Receptor_indexed(:,1);
y2 = Receptor_indexed(:,4)-Receptor_indexed(:,2);
norm_donor1 = [x1 y1 x2 y2];

d1 = Donor(:,3)-Donor(:,1);
d2 = Donor(:,4)-Donor(:,2);
r1 = Receptor(:,3)-Receptor(:,1);
r2 = Receptor(:,4)-Receptor(:,2);
norm_donor2 = [d1 d2 r1 r2];

z = zeros(length(x),1);
for j = 1:length(bacteria)
    x1 = find(strcmp(x,bacteria(j)));
    z(x1) = j;
end


%colormap
% x2 = linspace(0,953,4);
load('/Users/sophiesoham/Desktop/Eleni/ALL_ANALYSIS/Matlab variables/RWBmap2.mat');
y1 = 1:length(norm_donor1);
imagesc(recep);colormap(RWBmap);
caxis([0 7]);
hold on;
% 
% 
[coeff, score, latent, explained, tsquared] = princomp(alldata);
% Get the eigenvalues of B.
[A,L] = eig(coeff);
% Find the ordering largest to smallest.
[vals, inds] = sort(diag(L));
inds = flipud(inds);
vals = flipud(vals);
% Re-sort based on these.
A = A(:,inds);
L = diag(vals);
X = A(:,1:2)*diag(sqrt(vals(1:2)));
ind6 = find(coeff(1,:) <= 0.3);
ind9 = find(coeff(1,:) >= 0.5);
plot(X(ind6,1),X(ind6,2),'x',X(ind9,1),X(ind9,2),'o')
legend({'Topic 6';'Topic 9'})



figure;
scatter(coeff(:,1), coeff(:,2), 'filled', 'g'); hold on
scatter(coeff(:,2), coeff(:,3),'filled',   'r'); hold on
scatter(coeff(:,3), coeff(:,4),'filled', 'b'); hold on
hold off

a = rotatefactors(coeff(:,1:2),'Method','equamax');
figure;
biplot(a);


% %%
% 
% for i = 1: length(Index)
%     Don = Index(i).donor;
%     Rec = Index(i).receptor;
%     xs = find(Don(:,1) == 0);
%     xs2 = find(Don(:,3) == 0);
%     xs3 = find(Rec(:,4) <median(Rec(:,4)));
%     I = union(xs,xs2);
%     I1 = union(I, xs3);
%     Don(I1,:) = [];
%     Rec(I1,:) = [];
%        
%    Index(i).donor_corr1 = Don;
%    Index(i).receptor_corr1 = Rec;
%     
%         
% end
% 
% blah = [];
% blah1 = [];
% for i= 1:length(Index)
%     blah = [blah; Index(i).donor_corr1];
%     Donor_indexed = blah;
%     blah1 = [blah1; Index(i).receptor_corr1];
%     Receptor_indexed = blah1;
% end
% clear blah blah1
% 
% x1 = Donor_indexed(:,3)-Donor_indexed(:,1);
% y1 = Donor_indexed(:,4)-Donor_indexed(:,2);
% x2 = Receptor_indexed(:,3)-Receptor_indexed(:,1);
% y2 = Receptor_indexed(:,4)-Receptor_indexed(:,2);
% norm_donor3 = [x1 y1 x2 y2];
% 
% vari = var(norm_donor1');
% [v indx] = sort(vari', 'descend');
% 
% 
% x = [2 4 5 6 9 11];
% figure;
% for i = 1:length(x)
% subplot(3,2,i); bar(norm_donor1(x(i),:));ylim([-4 4])
% end
% 
% %Shannon Index
% data = all;
% ln_data = real(log(all));
% blah = -(data .* ln_data);
% for i = 1:size(blah,1)
%     for j = 1:size(blah,2)
%         if isnan(blah(i,j))== 1
%             blah(i,j) = 0;
%         end
%     end
% end
% alpha = 1./(blah);
% 

imagesc(x1);caxis([0 6]);colormap(RWBmap);
