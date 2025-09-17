using LinearAlgebra
using SparseArrays

header(f) = println(f,
"""
\\documentclass[dvisvgm,tikz]{standalone}
\\begin{document}
\\begin{tikzpicture}
""")

footer(f) = println(f,
"""
\\end{tikzpicture}
\\end{document}
""")

line(f, p1, p2) =
    println(f, "\\draw ($(p1[1]), $(p1[2])) -- ($(p2[1]),$(p2[2]));")

line(f, p1, p2, attr) =
    println(f, "\\draw[$attr] ($(p1[1]), $(p1[2])) -- ($(p2[1]),$(p2[2]));")

edge(f, n1, n2, attr) =
    println(f, "\\draw[$attr] ($n1) -- ($n2);")
    
node(f, name, p, attr, label) =
    println(f, "\\node[$attr] ($name) at ($(p[1]),$(p[2])) {$label};")

text(f, p, attr, label) =
    println(f, "\\node[$attr] at ($(p[1]),$(p[2])) {$label};")

# Set up picture coordinates
xy = [0.0 1.0 2.0 2.0 0.0 0.0 2.0 4.0 4.0 1.0 3.0;
      0.0 0.0 0.0 1.0 1.0 2.0 2.0 1.0 2.0 3.0 3.0]

# Set up sparse matrix
e = [1 2; 1 5; 2 3; 3 4; 4 5; 5 6; 5 7; 6 7; 6 10; 7 9; 7 10; 7 11; 8 9; 9 11; 10 11]
et = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 9; 8 9; 9 10; 10 11]
LA = sparse(e[:,2], e[:,1], ones(15), 11,11)
A = 0.25*(LA+LA') + Matrix{Float64}(I,11,11)
F = cholesky(A)
L = F.L

nstyle = "shape=rectangle, rounded corners, draw, align=center, top color=white, bottom color=blue!20"
nstyle2 = "shape=rectangle, rounded corners, draw, align=center, top color=white, bottom color=black!20, text=black!50"

open("sgraph-A.tex", "w") do f
    header(f)
    for j=1:11
        node(f, "n$j", xy[:,j], nstyle, j)
    end
    for i=1:15
        edge(f, "n$(e[i,1])", "n$(e[i,2])", "thick")
    end
    footer(f)
end

AA = copy(A)
for k = 1:10

    # Step of Cholesky
    AA[k,k] = sqrt(AA[k,k])
    AA[k+1:end,k] ./= AA[k,k]
    AA[k+1:end,k+1:end] .-= AA[k+1:end,k]*AA[k+1:end,k]'

    # Print current state of elimination
    open("sgraph-A-$k.tex", "w") do f
        header(f)

        # Nodes
        for j=1:k      node(f, "n$j", xy[:,j], nstyle2, j) end
        for j=k+1:11   node(f, "n$j", xy[:,j], nstyle, j)  end

        # L edges
        for j=1:11
            for i=j+1:11
                if j <= k
                    if A[i,j] != 0
                        edge(f, "n$j", "n$i", "thick,black!20,->")
                    elseif AA[i,j] != 0
                        edge(f, "n$j", "n$i", "thick,black!20,dashed,->")
                    end
                else
                    if A[i,j] != 0
                        edge(f, "n$j", "n$i", "ultra thick")
                    elseif AA[i,j] != 0
                        edge(f, "n$j", "n$i", "ultra thick,red,dashed")
                    end
                end
            end
        end
        
        footer(f)
    end

end

open("sgraph-L.tex", "w") do f
    header(f)
    for j=1:11
        node(f, "n$j", xy[:,j], nstyle, j)
    end
    for i=1:15
        edge(f, "n$(e[i,1])", "n$(e[i,2])", "thick,->")
    end
    for j=1:11
        for i = j+1:11
            if A[i,j] == 0 && L[i,j] != 0
                edge(f, "n$j", "n$i", "thick,red,dashed,->")
            end
        end
    end
    footer(f)
end

open("sgraph-LT.tex", "w") do f
    header(f)
    for j=1:11
        node(f, "n$j", xy[:,j], nstyle, j)
    end
    for i=1:15
        edge(f, "n$(e[i,1])", "n$(e[i,2])", "->")
    end
    for j=1:11
        for i = j+1:11
            if A[i,j] == 0 && L[i,j] != 0
                edge(f, "n$j", "n$i", "red,dashed,->")
            end
        end
    end
    for i = 1:10
        if A[et[i,2],et[i,1]] != 0
            edge(f, "n$(et[i,1])", "n$(et[i,2])", "ultra thick,->")
        else
            edge(f, "n$(et[i,1])", "n$(et[i,2])", "ultra thick,red,dashed,->")
        end
    end
    footer(f)
end

open("sgraph-AL-fig.tex", "w") do f
    println(f, "\$\$")
    println(f, "\\begin{bmatrix}")
    for i=1:11
        if A[i,1] != 0
            print(f, "\\times")
        end
        for j=2:11
            print(f, " & ")
            if A[i,j] != 0
                print(f, "\\times")
            end
        end
        if i != 11
            println(f, "\\\\")
        end
    end
    println(f, "\\end{bmatrix}, \\quad")
    println(f, "\\begin{bmatrix}")
    for i=1:11
        if A[i,1] != 0
            print(f, "\\times")
        elseif L[i,1] != 0
            print(f, "*")
        end
        for j=2:11
            print(f, " & ")
            if A[i,j] != 0
                print(f, "\\times")
            elseif L[i,j] != 0
                print(f, "*")
            end
        end
        if i != 11
            println(f, "\\\\")
        end
    end
    println(f, "\\end{bmatrix}")
    println(f, "\$\$")    
end
