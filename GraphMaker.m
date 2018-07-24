X = csvread("Results.csv", 2);
gens = X(:,1); best_host = X(:,2); bh_para = X(:,3); 
best_para = X(:,4); bp_host = X(:,5);

figure
plot(gens, best_host, gens, bh_para);
xlabel("Generations"); ylabel("Fitness");
legend("Best Host", "Best Host's Parasite");
set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);
