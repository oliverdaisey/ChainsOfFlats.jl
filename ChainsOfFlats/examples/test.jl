using .ChainsOfFlats
using Oscar

n=7
k=5
TropL = ChainsOfFlats.trivalent_tropical_linear_space(n,k, true)
bases = [I for (I,pI) in zip(pluecker_indices(TropL), tropical_pluecker_vector(TropL)) if !iszero(pI)]

E = identity_matrix(QQ, n)
matroidPolytope = convex_hull([E[first(b),:] + E[last(b),:] for b in bases])
matroidVertices = vertices(matroidPolytope)
matroidFacets = IncidenceMatrix(faces(matroidPolytope, dim(matroidPolytope)-1))
matroidFacets = [findall(matroidFacets[i,:]) for i in 1:nrows(matroidFacets)]
for matroidFacet in matroidFacets
    println(sum([ matroidVertices[i] for i in matroidFacet]))
end