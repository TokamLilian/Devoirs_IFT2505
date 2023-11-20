mutable struct SimplexTableau
    z_c     ::Array{Float64}  # coûts réduits
    Y       ::Array{Vector{Float64}}  # contraintes linéaires
    x_B     ::Array{Float64}  # solution de base
    obj     ::Float64         # valeur de la fonction objectif
    b_idx   ::Array{Int64}    # indices des variables de bases
  
    function SimplexTableau(z_c, Y, x_B, obj, b_idx)
      new(z_c, Y, x_B, obj, b_idx)
    end
  
  end 
  
  
  function draw_table(t::SimplexTableau; verbose::Bool)
    println("------------------------------")
    A = hcat(t.Y...)
    M = [A' t.x_B; t.z_c' t.obj]
    verbose && display(M)
    println("------------------------------")
  end
  
  
  function pivot_point(tableau_simplexe)
    # Récupération des indices des variables de base
    b_idx = tableau_simplexe.b_idx
  
    # Calcul des coûts réduits des variables non de base
    z_c = tableau_simplexe.z_c
  
    #Recuperer le nombre de variables de départ
    nombre_variables = length(z_c) - length(b_idx)
  
    # Sélection de la variable d'entrée (variable avec le coût réduit le plus négatif)
    if 1 in b_idx
      indice_variable_entree = argmax(z_c[1:nombre_variables])
    else
      indice_variable_entree = argmin(z_c)
    end
  
    # Récupération des coefficients de la colonne correspondant à la variable d'entrée
    colonne_variable_entree = [ligne[indice_variable_entree] for ligne in tableau_simplexe.Y]
  
    # Sélection de la variable de sortie (variable de base avec le rapport positif minimal)
    indices_sortants = findall(x -> x != 0, colonne_variable_entree)    ## specifique au probleme du devoir 3
  
    rapport_minimal = Inf
    indice_variable_sortante = nothing
  
    for indice_sortant in indices_sortants
        rapport = tableau_simplexe.x_B[indice_sortant] / colonne_variable_entree[indice_sortant]
        if rapport < rapport_minimal
            rapport_minimal = rapport
            indice_variable_sortante = indice_sortant
        end
    end
  
    # Retourner les indices de la variable d'entrée et de sortie
    return indice_variable_entree, b_idx[indice_variable_sortante]
  end
  
  
  function pivoting!(t::SimplexTableau)
      m, n = length(t.Y), length(t.Y[1])
  
      entering, exiting = pivot_point(t)
      exiting_temp = exiting
  
      #Recuperer la position(ligne) de la variable de sortie
      exiting = findfirst(x -> x==exiting, t.b_idx)   
  
      @show (entering, exiting)
  
      # Pivoting: exiting-row, entering-column
      # updating exiting-row
      coef = t.Y[exiting][entering]
      t.Y[exiting, :] /= coef
      t.x_B[exiting] /= coef
  
      # updating other rows of Y
      for i in setdiff(1:m, exiting)
        coef = t.Y[i][entering]
        t.Y[i, :] -= coef * t.Y[exiting, :]
        t.x_B[i] -= coef * t.x_B[exiting]
      end
  
      # updating the row for the reduced costs
      coef = t.z_c[entering]
      t.z_c[:] -= coef * t.Y[exiting][:]
      t.obj -= coef * t.x_B[exiting]
  
      # Updating b_idx
      t.b_idx[findall(x -> x == exiting_temp, t.b_idx)] .= entering
      return t
  
  end 
  
  #feasible prend un tableau du simplex et retourne si il contient une solution primal rélisable
  function feasible(t::SimplexTableau)
    # Vérifier les couts reduits
    couts_reduits = all(t.z_c[:, end] .>= 0)
  
    # Vérifier la non-négativité des variables
  
    non_negativite_variables = all(t.x_B[:, 1:end] .>= 0)
  
    # La solution est réalisable si les deux conditions sont satisfaites
    return non_negativite_variables
  end

  using LinearAlgebra
nouveauTableau = SimplexTableau([0, 3, 4, 5, 0, 0, 0], 
                          [ [-1//1, 0, -1, 0, 1, 0, 0],
                            [0, -1, -2, -1, 0, 1, 0],
                            [0,  0,  0, -1, 0, 0, 1] ],
                          [-3, -4, -2],
                          0,
                          [5, 6, 7])

function DualsimplexSolver(t::SimplexTableau; verbose::Bool = false)
  k = 1
  nmax = 1000

  while !feasible(t) && k < nmax
      
      k+=1
      draw_table(t, verbose=true)
      pivoting!(t)
    end
    
  draw_table(t, verbose=true)
  return t.obj
end                    

DualsimplexSolver(nouveauTableau, verbose = true)

