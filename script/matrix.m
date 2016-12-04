matrix[l_, \[Alpha]_, \[Beta]_, \[Kappa]_, \[Tau]_, t1_,t2_] := 
    Normal[MatrixExp[(t1+I*t2)*N[getmatrix[l, \[Alpha], \[Beta], \[Kappa], \[Tau]],100]]]
 
getmatrix[l_, \[Alpha]_, \[Beta]_, \[Kappa]_, \[Tau]_] := 
    (-\[Kappa])*getmatrixE2[l] - \[Tau]*getmatrixE3[l] + 
     ((-\[Alpha]^(-1) + 1/\[Beta])/(2*\[Alpha]))*getmatrixE3[l] . 
       getmatrixE3[l] + (1/(2*\[Alpha]))*getmatrixLaplace[l]
 
getmatrixE2[l_] := (1/2)*SparseArray[If[l > 0, 
       {Band[{1, 2}] -> Table[-Sqrt[(l - j)*(l + j + 1)], {j, -l, l - 1}], 
        Band[{2, 1}] -> Table[Sqrt[(l - j)*(l + j + 1)], {j, -l, l - 1}]}, 
       {1, 1} -> 0], {2*l + 1, 2*l + 1}]
 
getmatrixE3[l_] := I*SparseArray[Band[{1, 1}] -> Table[j, {j, -l, l}], 
      {2*l + 1, 2*l + 1}]
 
getmatrixLaplace[l_] := (-l)*(l + 1)*SparseArray[Band[{1, 1}] -> 1, 
      {2*l + 1, 2*l + 1}]

