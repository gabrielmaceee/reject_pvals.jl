using MultipleTesting
using Plots

# ajouter un seuil, et retourner les (indexs) pvaleurs corrigées inférieures à ce seuil.
# changer l'appel : str au lieu de fonction
# ajouter un graph avec les p-vals non corrigées et le seuil corrigées
# sur le graph : différencier les significatifs et les non 

# A partir d'un vecteur de p-valeurs, d'une fonction méthode de MultipleTesting, et d'un seuil, 
# retourne l'index des p-valeurs inférieures au seuil, et les p-valeurs correspondantes.
function adjust_select(pvals, methode::String, seuil::Number)
    pvals_adj = p_adjust(pvals, methode)
    index = findall(pvals_adj .< seuil)
    return index,  pvals_adj[index]
end

function adjust_select(pvals, methode::String, seuil::Vector{<:Number})
    pvals_adj = p_adjust(pvals, methode)
    index = Vector()
    for i in 1:size(pvals)[1]
        if pvals_adj[i] < seuil[i] 
            index = [index; i]
        end
    end
    return index,  pvals_adj[index]
end

function adjust_select(pvals, methode, seuil::Number)
    pvals_adj = p_adjust(pvals, to_str(methode))
    index = findall(pvals_adj .< seuil)
    return index,  pvals_adj[index]
end

function adjust_select(pvals, methode, seuil::Vector{<:Number})
    pvals_adj = p_adjust(pvals, to_str(methode))
    index = Vector()
    for i in 1:size(pvals)[1]
        if pvals_adj[i] < seuil[i] 
            index = [index; i]
        end
    end
    return index,  pvals_adj[index]
end

# Applique la correction de Simes à un vecteur de p-valeurs, retourne les p-valeurs corrigées.
function simes(pvals::Vector{<:Number})
    m = size(pvals)[1]
    i = 1:m
    correction =  m ./ i
    pvals_adj = correction .* sort(pvals)
    mat = hcat(pvals, 1:m)
    mat = sortslices(mat, dims=1, by = x->x[1])   
    return pvals_adj[Int.(mat[:,2])]
end

# Retourne le nom de la méthode de MultipleTesting
function to_str(methode)
    if methode == Bonferroni()
        return "Bonferroni"
    elseif  methode == Hochberg()
        return "Hochberg"
    elseif methode == Holm()
        return "Holm"
    else 
        return "Nop"
    end
end


# Retourne le nom de la méthode de MultipleTesting
function p_adjust(pvals::Vector{<:Number}, methode::String)
    if methode == "Bonferroni"
        return adjust(pvals, Bonferroni())
    elseif  methode == "Hochberg"
        return adjust(pvals, Hochberg())
    elseif methode == "Holm"
        return adjust(pvals, Holm())
    elseif methode == "Simes"
        return simes(pvals)
    else 
        return "Nop"
    end
end

# A partir d'un nom de méthode, d'une taille de vecteurs de p-valeus et d'un seuil, retourne le seuil corrigés
function get_seuil(methode::String, m::Int, risk::Number)
    if methode == "Bonferroni"
        return repeat([risk / m],m)
    elseif  methode == "Simes"
        seuil = Vector()
        for i in 1:m
            seuil = [seuil; (risk * i) ./ m]     
        end
        return seuil
    elseif  methode == "Hochberg"
        i = 1:m 
        return risk ./ (m .- i.+ 1)
    elseif methode == "Holm"
        i = 1:m 
        return risk ./ (m .- i .+ 1)
    else
        return "Nop"
    end
end

# A partir d'un vecteur de p-valeurs et d'un seuil, retourne si les p-valeurs sont inférieurs au seuil.
function get_rejet(pvals::Vector{<:Number}, seuil, method::String, aj = Bool)  
    if aj == true
        rejet = adjust_select(pvals, method, seuil)[1]
    else
        rejet = Vector()
        for i in 1:size(pvals)[1]
            if pvals[i] < seuil[i] 
                rejet = [rejet; i]
            end
        end
    end

    r = Vector()
    for i in 1:size(pvals)[1]
        if i in rejet
            r = [r; true]
        else r = [r; false]
        end
    end
    return r
end




# Graph où les seuils sont corrigés, mais pas les p-valeurs :
# multiple dispatching en fonction du type de seuil, si seuil = une seule valeur : seuil à corriger :
function graph_pvals(pvals::Vector{<:Number},  seuil::Number, method::String)
    pvals = sort(pvals)
    seuils = get_seuil(method, size(pvals)[1], seuil)
    colors = Int.(get_rejet(pvals, seuils, method, false))
    x = 1:size(pvals)[1]
    pvals_accept = pvals[colors .== 1]
    x_accept = x[colors .== 1]

    pvals_reject = pvals[colors .== 0]
    x_reject = x[colors .== 0]
    scatter(x_accept, pvals_accept, color=:blue, label="Acceptés (p > seuil)", legend=:top)
    scatter!(x_reject, pvals_reject, color=:red, label="Rejetés (p <= seuil)")
    #scatter(1:size(pvals)[1], sort(pvals), label = "p-values", color = Int.(get_rejet(pvals, seuil)))
    plot!(x, seuils, label = "Seuil")
end

function graph_pvals(pvals::Vector{<:Number}, seuil::Vector{<:Number}, method::String)
    colors = Int.(get_rejet(pvals, seuil, method, false))
    x = 1:size(pvals)[1]
    pvals_accept = pvals[colors .== 1]
    x_accept = x[colors .== 1]

    pvals_reject = pvals[colors .== 0]
    x_reject = x[colors .== 0]
    scatter(x_accept, pvals_accept, color=:blue, label="Acceptés (p > seuil)", legend=:top)
    scatter!(x_reject, pvals_reject, color=:red, label="Rejetés (p <= seuil)")
    #scatter(1:size(pvals)[1], sort(pvals), label = "p-values", color = Int.(get_rejet(pvals, seuil))) # , ylim=(0, 0.05)
    plot!(x, seuil, label = "Seuil")
end

# Graph où les p-valeurs sont corrigées :
function graph_pvals_corr(pvals::Vector{<:Number},  risk::Number, method::String)
    seuils = repeat([risk],size(pvals)[1])
    colors = Int.(get_rejet(pvals, seuils, method, true))
    pvals = p_adjust(pvals, method)
    x = 1:size(pvals)[1]
    pvals_accept = pvals[colors .== 1]
    x_accept = x[colors .== 1]

    pvals_reject = pvals[colors .== 0]
    x_reject = x[colors .== 0]
    scatter(x_accept, pvals_accept, color=:blue, label="Acceptés (p > seuil)", legend=:top)
    scatter!(x_reject, pvals_reject, color=:red, label="Rejetés (p <= seuil)")
    plot!(x, seuils, label="Seuil")
end



