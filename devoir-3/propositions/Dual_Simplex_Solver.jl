function pivot!(M::Matrix, i::Int, j::Int)
    m, n = size(M)
    @assert M[i, j] != 0
    M[i, :] = M[i, :]/M[i, j]
    for k in setdiff(1:m, i)
        M[k, :] -= M[k, j] * M[i, :]
    end
    return M
end

function getxB(M::Matrix)
    m, n = size(M)
    return M[1:m-1, end]
end

function isOneHot(v::Vector)
    n = length(v)
    return (sum(iszero, v) == n-1) && (sum(isone, v) == 1)
end 

function findInitialBasis!(M::Matrix)
    m, n = size(M)
    m-=1
    n-=1
    basis = [-1 for _ in 1:m]
    for i in 1:n
        if isOneHot(M[1:end-1, i])
            index = findfirst(isone, M[:, i])
            basis[index] = i
        end
    end
    @assert !any(t-> t == -1, basis) "problem not caconical"
    for i in 1:m
        j = basis[i]
        M = pivot!(M, i, j)
    end
    return basis
end

#feasible prend un tableau du simplex et retourne si il contient une solution primal rélisable
function feasible(M::Matrix{T}) where T
    # Vérifier les contraintes linéaires
    contraintes_lineaires = all(M[:, end] .>= 0)

    # Vérifier la non-négativité des variables
    non_negativite_variables = all(M[:, 1:end-1] .>= 0)

    # La solution est réalisable si les deux conditions sont satisfaites
    return contraintes_lineaires && non_negativite_variables
end

#retourrne la ligne du tableau du simplex dont la variable de base va sortir de la base
function findLeavingVarInDual(M::Matrix{T}) where T
    variables_sortie_positives = findall(x -> x > 0, M[:, end-1])

    # S'il n'y a pas de variables de sortie positives, le problème n'est pas borné
    isempty(variables_sortie_positives) && throw(ArgumentError("Problème non borné"))

    ratios = [(i, M[i, end] / M[i, argmax(M[i, 1:end-2])]) for i in variables_sortie_positives]
    
    # Trouver la ligne correspondant au ratio minimum
    ligne_sortie = argmin([x[2] for x in ratios])

    return ligne_sortie
end

#retourne la variable qui va entrer dans la base. Si aucune variable ne peut entrer dans
#la base, on retourne -1.
function findEnteringVarInDual(M::Matrix{T}, leaving::Int) where T
    coefficients_premiere_ligne = M[1, 1:end-1]

    # Vérifier si tous les coefficients de la première ligne sont négatifs ou nuls
    if all(x -> x <= 0, coefficients_premiere_ligne)
        return -1  # Aucune variable ne peut entrer dans la base
    else
        # Trouver l'indice de la variable avec le coefficient le plus négatif
        indice_variable_entree = argmin(coefficients_premiere_ligne)
        return indice_variable_entree
    end
end


function DualsimplexSolver(A::Matrix{T}, b::Vector, c::Vector; verbose::Bool = false) where T
    @assert all(c .>= 0) #dual feasibility
    M = [A b; c' 0]
    basis = findInitialBasis!(M)
    k = 1
    nmax = 1000
    while !feasible(M) && k < nmax
        k+=1
        verbose && display(M)
        leaving = findLeavingVarInDual(M)
        entering = findEnteringVarInDual(M, leaving)
        if entering == -1
            println("aucune variable ne peut entrer dans la base.")
            println("problème non réalisable")
        end
        verbose && @show (entering, leaving)
        basis[leaving] = entering
        M = pivot!(M, leaving, entering)
    end
    verbose && display(M)
    m, n = size(M)
    xstar = zeros(T, n - 1)
    xstar[basis] = getxB(M)
    return xstar
end


using LinearAlgebra

Atilde = [1//1 0 1 0; 0 1 2 1; 0 0 0 1]
A = -[Atilde -I(3)]
@show size(A)
b = -[3, 4, 2]
@show size(b)
c = [0, 3, 4, 5, 0, 0, 0]
@show size(c);

xstar = DualsimplexSolver(A, b, c, verbose=true);
#xstar = simplexSolver(A, b, c, verbose=true);