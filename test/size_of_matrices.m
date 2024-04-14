%% increase k
figure(1)
clf
tl = tiledlayout(2, 2);
nexttile;
hold on
n = 2;
m = 2;
T = 12;

klist = 1:6;
ss = [];
ds = [];
fs = [];
c = linspecer(3);
for i = 1:length(klist)
    ss(i) = sparse_size(n, m, T, i);
    ds(i) = dense_size(n, m, T, i);
    fs(i) = full_size(n, m, T, i);
end

plot(klist, fs, '-o', 'color', c(1, :));
plot(klist, ds, '-o', 'color', c(2, :));
plot(klist, ss, '-o', 'color', c(3, :));
legend({'Full', 'Dense', 'Sparse'}, 'location', 'northwest')
xlabel('k (degree)')
ylabel('PSD size')
set(gca, 'YScale', 'log')
title(sprintf('PSD Size vs. Degree (n=%d, m=%d, T=%d)', n, m, T), 'fontsize', 14)
grid on

%% increase T
nexttile
hold on
n = 2;
m = 2;
k = 2;

Tlist = 4:20;
ss = [];
ds = [];
fs = [];
c = linspecer(3);
for i = 1:length(Tlist)
    ss(i) = sparse_size(n, m, i, k);
    ds(i) = dense_size(n, m, i, k);
    fs(i) = full_size(n, m, i, k);
end

plot(Tlist, fs, '-o', 'color', c(1, :));
plot(Tlist, ds, '-o', 'color', c(2, :));
plot(Tlist, ss, '-o', 'color', c(3, :));
legend({'Full', 'Dense', 'Sparse'}, 'location', 'northwest')
xlabel('T (Number of Samples)')
ylabel('PSD size')
set(gca, 'YScale', 'log')
title(sprintf('PSD Size vs. Number of Samples (n=%d, m=%d, k=%d)', n, m, k), 'fontsize', 14)
grid on

%% increase n
nexttile
hold on
T = 12;
m = 2;
k = 2;

nlist = 2:8;
ss = [];
ds = [];
fs = [];
c = linspecer(3);
for i = 1:length(nlist)
    ss(i) = sparse_size(i, m, T, k);
    ds(i) = dense_size(i, m, T, k);
    fs(i) = full_size(i, m, T, k);
end

plot(nlist, fs, '-o', 'color', c(1, :));
plot(nlist, ds, '-o', 'color', c(2, :));
plot(nlist, ss, '-o', 'color', c(3, :));
legend({'Full', 'Dense', 'Sparse'}, 'location', 'northwest')
xlabel('n (Number of States)')
ylabel('PSD size')
set(gca, 'YScale', 'log')
title(sprintf('PSD Size vs. Number of States (m=%d, k=%d, T=%d)', m, k, T), 'fontsize', 14)
grid on

%% increase m
nexttile
hold on
T = 12;
m = 2;
k = 2;

mlist = 1:8;
ss = [];
ds = [];
fs = [];
c = linspecer(3);
for i = 1:length(mlist)
    ss(i) = sparse_size(n, i, T, k);
    ds(i) = dense_size(n, i, T, k);
    fs(i) = full_size(n, i, T, k);
end

plot(mlist, fs, '-o', 'color', c(1, :));
plot(mlist, ds, '-o', 'color', c(2, :));
plot(mlist, ss, '-o', 'color', c(3, :));
legend({'Full', 'Dense', 'Sparse'}, 'location', 'northwest')
xlabel('m (Number of Inputs)')
ylabel('PSD size')
set(gca, 'YScale', 'log')
title(sprintf('PSD Size vs. Number of Inputs (n=%d, k=%d, T=%d)', n, k, T), 'fontsize', 14)
grid on

function z = sparse_size(n, m, T, k)
    nvar = (n+m);
    z = nchoosek(nvar+k, k);
end

function z = dense_size(n, m, T, k)
    nvar = n*(n+m);
    z = nchoosek(nvar+k, k);
end

function z = full_size(n, m, T, k)
    nvar = n*(n+m) + n*T + m*(T-1);
    z = nchoosek(nvar+k, k);
end