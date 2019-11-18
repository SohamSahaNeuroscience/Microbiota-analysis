clearvars -except x

load('alldata.mat');% copied from the excel sheet
load('allotuids.mat'); % ["Phylum","Class","Order","Family"]


colid = alldata(end,:);
[i, m] = sort(colid, 'ascend');

data = alldata(:,m);

ref = (21:62);

finaldata = data(1:end-1,ref);
finaldata(:,19) = [];% values are too low
q = sum(finaldata,2);
z = find(q == 0);%find all zeo values
finaldata(z,:) = [];
clear q z

a = finaldata;
%min-max normalization
for i = 1:size(a,1)
    b = a(i,:);
    finaldata_n(i,:) = (a(i,:)-min(b))/(max(b)-min(b));
    clear b
end

s = sum(finaldata_n,2);
e = find(isnan(s)==1);
finaldata_n(e,:) = [];
clear s e

b = cellstr(allotuids);

b1 = unique(b(1:end-1,4));%family
b2 = unique(b(1:end-1,1));%phylum
b3 = unique(b(1:end-1,6));%genus
%FOR FAMILY
for i = 1:length(b1)
Index(i).family = find(strcmp(b(1:length(finaldata_n),5), b1(i)));
Index(i).fam = finaldata_n(Index(i).family,:);
end
%copy b1 or b2 physically into Index files

%FOR PHYLUM
for i = 1:length(b2)
Index(i).phylum = find(strcmp(b(1:length(finaldata_n),2), b2(i)));
Index(i).phy = finaldata_n(Index(i).phylum,:);
end

%FOR GENUS
for i = 1:length(b3)
Index(i).genus = find(strcmp(b(1:length(finaldata_n),6), b3(i)));
Index(i).gen = finaldata_n(Index(i).family,:);
end


%MEAN for every family or phylum
for i = 1:size(Index,2)
    bla = Index(i).fam;
    bla1 = Index(i).phy;
    bla2 = Index(i).gen;
    
    m = nanmean(bla,1);
    m1= nanmean(bla1,1);
    m2 = nanmean(bla2,1);
    
    Index(i).fam_mean = m;
    Index(i).phy_mean = m1;
    Index(i).gen_mean = m2;
    clear bla bla1 m m1 bla2 m2
end
clear i

c = [];
c1 = [];
c2 = [];
for i = 1:size(Index,2)
c = [c; Index(i).fam_mean];
c1 = [c1; Index(i).phy_mean];
c2 = [c2; Index(i).gen_mean];
end

l = find(isnan(sum(c,2)) == 1);
l1 = find(isnan(sum(c1,2)) == 1);
l2 = find(isnan(sum(c2,2)) == 1);

c(l,:) = [];
c1(l1,:) = [];
c2(l2,:) = [];


c = Index(1).FAMILY% = c;
c1 = Index(1).PHYLUM% = c1;
Index(1).GENUS = c2;


%Shannon Index
[H,VarH]=index_SaW(c,2);
[H1,VarH1]=index_SaW(c1,2);
[H2,VarH2]=index_SaW(c2,2);

%Sort data

v = mean(c');
[vs vi] = sort(v, 'descend');
c_sort = c(vi,:);
b1_sort = b1(vi,:);


v1 = mean(c1');
[vs1 vi1] = sort(v1, 'descend');
c1_sort = c1(vi1,:);
b2_sort = b2(vi1,:);

v2 = mean(c2');
[vs2 vi2] = sort(v2, 'descend');
c2_sort = c2(vi2,:);
b3_sort = b3(vi2,:);

a = c_sort(1,:);%family
aa = c1_sort(1,:);%phylum
aaa = c2_sort(1,:);%genus
figure; scatter(a,aa); hold on; scatter(aa, aaa); hold off


%autocorrelation
p = Index(1).PHYLUM;
p(:,19:31) = [];
p(:,1:8) = [];
cr_phy = corr(p, p);
%phyu = triu(cr_phy);
f = Index(1).FAMILY;
f(:,19:31) = [];
f(:,1:8) = [];
cr_fam = corr(f, f);
%phyl = tril(cr_fam);
cr_gen = corr(Index(1).GENUS, Index(1).GENUS);
%phyl = tril(cr_gen);


ab_fam = Index(1).FAMILY(:,19:31);
cr_abx_fam = corr(ab_fam, ab_fam);
ab_phy = Index(1).PHYLUM(:,19:31);
cr_abx_phy = corr(ab_phy, ab_phy);

figure;imagesc(cr_fam);colormap(RWBmap);caxis([0 0.65])
figure;imagesc(cr_phy);colormap(RWBmap);caxis([0 0.65])
figure;imagesc(cr_abx_fam);colormap(RWBmap);caxis([0 1])
figure;imagesc(cr_abx_phy);colormap(RWBmap);caxis([0 1])

mat = phyu+phyl;

for i = 1:size(mat,1)
    for j = 1:size(mat,2)        
        if mat(i,j) < 0.4
            mat(i,j) = 0;
        end
    end
end
