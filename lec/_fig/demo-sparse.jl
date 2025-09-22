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

# Set up graph
xy = [3.0 0.0 1.0 2.0 2.0 0.0 4.0 5.0 4.0 5.0 0.0 2.0 4.0 2.0 4.0;
      4.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0 1.0 2.0 2.0 2.0 3.0 3.0]
e = [1  1  2 2 3 4 5 5  6  7 7 8  9  9  10 11 11 12 12 13;
     14 15 3 6 4 5 6 12 11 8 9 10 10 13 13 12 14 13 14 15]

n = size(xy)[2]
m = size(e)[2]

# Set up sparse matrix and factor
LA = sparse(e[2,:], e[1,:], ones(m), n,n)
A = 0.25*(LA+LA') + Matrix{Float64}(I,n,n)
F = cholesky(A)
L = F.L

# Build etree
et = zeros(Int,2,n)
for j=1:n-1
    et[1,j] = j
    et[2,j] = j+findfirst(L[j+1:end,j] .!= 0)
end

nstyle = "shape=rectangle, rounded corners, draw, align=center, top color=white, bottom color=blue!20"
nstyle2 = "shape=rectangle, rounded corners, draw, align=center, top color=white, bottom color=black!20, text=black!50"

nclosed(c) = "\\node{{\\color{$c}\$\\bullet\$}};"
nopen(c) = "\\node{{\\color{$c}\$\\circ\$}};"
        
function Aspy(f, A)
    println(f, "\\matrix[draw,anchor=south west,column sep=-4pt,row sep=-4pt] at (6.0,-0.2) {")
    for i=1:n
        print(f, A[i,1])
        for j=2:n
            print(f, " & ")
            print(f, A[i,j])
        end
        println(f, "\\\\");
    end
    println(f, "};")
end

open("sgraph-A.tex", "w") do f
    header(f)
    for j=1:n
        node(f, "n$j", xy[:,j], nstyle, j)
    end
    for i=1:m
        edge(f, "n$(e[1,i])", "n$(e[2,i])", "ultra thick")
    end
    Aspy(f, [A[i,j] != 0 ? nclosed("black") : "" for i=1:n, j=1:n])
    footer(f)
end

AA = copy(A)
for k = 1:n-1

    # Step of Cholesky
    AA[k,k] = sqrt(AA[k,k])
    AA[k+1:end,k] ./= AA[k,k]
    AA[k+1:end,k+1:end] .-= AA[k+1:end,k]*AA[k+1:end,k]'

    # Print current state of elimination
    open("sgraph-A-$k.tex", "w") do f
        header(f)

        # Nodes
        for j=1:k     node(f, "n$j", xy[:,j], nstyle2, j) end
        for j=k+1:n   node(f, "n$j", xy[:,j], nstyle, j)  end

        # L edges
        for j=1:n
            for i=j+1:n
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

        function c(i,j)
            gscale = (i <= k || j <= k) ? "!40" : ""
            if A[i,j] != 0
                nclosed("black" * gscale)
            elseif AA[i,j] != 0
                nopen("red" * gscale)
            else
                ""
            end
        end
        Aspy(f, [c(i,j) for i=1:n, j=1:n])
        
        footer(f)
    end

end

open("sgraph-L.tex", "w") do f
    header(f)
    for j=1:n
        node(f, "n$j", xy[:,j], nstyle, j)
    end
    for i=1:m
        edge(f, "n$(e[1,i])", "n$(e[2,i])", "thick,->")
    end
    for j=1:n
        for i = j+1:n
            if A[i,j] == 0 && L[i,j] != 0
                edge(f, "n$j", "n$i", "thick,red,dashed,->")
            end
        end
    end
    
    function c(i,j)
        if A[i,j] != 0
            nclosed("black")
        elseif L[i,j] != 0 || L[j,i] != 0
            nopen("red")
        else
            ""
        end
    end
    Aspy(f, [c(i,j) for i=1:n, j=1:n])
    
    footer(f)
end

open("sgraph-LT.tex", "w") do f
    header(f)
    for j=1:n
        node(f, "n$j", xy[:,j], nstyle, j)
    end
    for i=1:m
        edge(f, "n$(e[1,i])", "n$(e[2,i])", "->")
    end
    for j=1:n
        for i = j+1:n
            if A[i,j] == 0 && L[i,j] != 0
                edge(f, "n$j", "n$i", "red,dashed,->")
            end
        end
    end
    for i = 1:n-1
        if A[et[2,i],et[1,i]] != 0
            edge(f, "n$(et[1,i])", "n$(et[2,i])", "ultra thick,->")
        else
            edge(f, "n$(et[1,i])", "n$(et[2,i])", "ultra thick,red,dashed,->")
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
        for j=2:n
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
    for i=1:n
        if A[i,1] != 0
            print(f, "\\times")
        elseif L[i,1] != 0
            print(f, "*")
        end
        for j=2:n
            print(f, " & ")
            if A[i,j] != 0
                print(f, "\\times")
            elseif L[i,j] != 0
                print(f, "*")
            end
        end
        if i != n
            println(f, "\\\\")
        end
    end
    println(f, "\\end{bmatrix}")
    println(f, "\$\$")    
end
