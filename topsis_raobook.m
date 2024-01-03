clc
clear
dmm = zeros (7:7);
m = size(dmm);
normal = zeros (7:7);
square_sum = zeros(7:7);
dmm(1,1) = 44; dmm(1,2) = 2.02; dmm(1,3) = 2.15; dmm(1,4) = 212; dmm(1,5) = 0.22; dmm(1,6) = 1007; dmm(1,7) = 276;
dmm(2,1) = 54; dmm(2,2) = 1.76; dmm(2,3) = 2.27; dmm(2,4) = 157; dmm(2,5) = 0.29; dmm(2,6) = 940; dmm(2,7) = 345;
dmm(3,1) = 61; dmm(3,2) = 1.69; dmm(3,3) =2.2; dmm(3,4) = 222; dmm(3,5) = 0.21; dmm(3,6) = 989; dmm(3,7) = 354;
dmm(4,1) = 90; dmm(4,2) = 3; dmm(4,3) = 2; dmm(4,4) = 250; dmm(4,5) = 0.2; dmm(4,6) = 950; dmm(4,7) = 500;
dmm(5,1) = 32; dmm(5,2) = 1.93; dmm(5,3) = 2.8; dmm(5,4) = 180; dmm(5,5) = 0.56; dmm(5,6) = 1485; dmm(5,7) = 48;
dmm(6,1) = 46; dmm(6,2) = 1.46; dmm(6,3) = 2.39; dmm(6,4) = 210; dmm(6,5) = 0.76; dmm(6,6) = 1666; dmm(6,7) = 199;
dmm(7,1) = 58; dmm(7,2) = 1.68; dmm(7,3) = 2.37; dmm(7,4) = 266; dmm(7,5) = 0.43; dmm(7,6) = 1450; dmm(7,7) = 233;
pij = dmm(7:7)./sum(dmm);
for j=1:7
    for i=1:7
        dmm_square(i,j) = (dmm(i,j).^2);
    end
end
cumsum = sum(sum(dmm_square));
weight = [0.28 0.14 0.05 0.24 0.19 0.05 0.05];
for j=1:7
    for i=1:7
        normal (i,j) = dmm(i,j)./sqrt(cumsum);
        weighted_normal (i,j) = normal(i,j).*weight(1,j);
    end
end

lnm = -1/log(size(dmm,1));
lnNormDmm = log(pij);
E = lnm.* sum(pij.*lnNormDmm);
dj = ones(1,size(E,2))-E;
weightEntropy = dj./sum(dj);
N = dmm(7:7)./repmat(dmm_square, [size(dmm,1)]);
Apositive = max(N);
Anegative = min(N);
ApositivMtrix = repmat(Apositive,size(N,1),1);
AnegativMtrix = repmat(Anegative,size(N,1),1);
s1 = (N-ApositivMtrix).^2;
s2 = (N-AnegativMtrix).^2;
for (j=1:1:size(s1,1))
    sumApositive(j) = sum(s1(j,:));
end
for (j=1:1:size(s2,1))
    sumAnegative(j) = sum(s2(j,:));
end
dpositive = sqrt(sumApositive);
dnegative = sqrt(sumAnegative);
sumD = dnegative + dpositive;
cc = dnegative./sumD;
y = sort(cc, 'descend');




