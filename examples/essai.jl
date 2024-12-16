using reject_pvals
using MultipleTesting
using Plots


pvaleurs = [0.01, 0.0008, 0.9, 0.7, 0.009, 0.05, 0.001, 0.0000000001, 0, 1]

adjust_select(pvaleurs, Bonferroni(), 0.05)

seuils = [0.1:0.1:1;]


simes(pvaleurs)

pvaleurs = [0.01, 0.0008, 0.09, 0.07, 0.009, 0.05, 0.001, 0.0000000001, 0, 0.2]

seuils = get_seuil("Hochberg", 10, 0.05)

get_seuil("Simes", size(pvaleurs)[1], 0.05)

graph_pvals_corr(pvaleurs, 0.05, "Hochberg")
graph_pvals_corr(pvaleurs, 0.05, "Holm")
graph_pvals_corr(pvaleurs, 0.05, "Bonferroni")
graph_pvals_corr(pvaleurs, 0.05, "Simes")



graph_pvals(pvaleurs,  0.05, "Hochberg")
graph_pvals(pvaleurs,  0.05, "Holm")
graph_pvals(pvaleurs,  0.05, "Bonferroni")
graph_pvals(pvaleurs,  0.05, "Simes")



p_adjust(pvaleurs, "Simes")





