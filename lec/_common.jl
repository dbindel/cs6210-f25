"""
    slide(title) do block end

Print a slide, contents printed by `block()`
"""
function slide(block, title)
    println("## $title\n")
    block()
    println("")
end


"""
    slatex() do block end

Print display LaTeX section, contents printed by `block()`
"""
function slatex(block)
    println("\$\$")
    block()
    println("\$\$")
end


"""
    svector(x) do i,xi f end

Print LaTeX vector, entry i given by `f(i,xi)`
"""
function svector(transform, x)
    n = length(x)
    println("\\begin{bmatrix}")
    for j = 1:n
        print(transform(j,x[j]))
        println(j < n ? " \\\\" : "")
    end
    println("\\end{bmatrix}")
end


"""
    smatrix(A) do i,j,xi aij end

Print LaTeX matrix, entry i,j given by `f(i,j,aij)`
"""
function smatrix(transform, A)
    m, n = size(A)
    println("\\begin{bmatrix}")
    for i = 1:m
        for j = 1:n
            print(transform(i,j,A[i,j]))
            print(j < n ? " & " : " ")
        end
        println(i < m ? " \\\\ " : "")
    end
    println("\\end{bmatrix}")
end


"""
    scolor(c,x)

Return a LaTeX string coloring text `x` with color `c`
"""
scolor(c, x) = "\\color{$c}{$x}"


"""
    print_spause()

Print a Quarto slide pause line
"""
print_spause() = println("\n. . .\n")


"""
    print_mini_spy(A; nzmark="\\times")

Print a LaTeX matrix spy plot, marking nonzeros with `nzmark`
"""
print_mini_spy(A; nzmark="\\times") =
    smatrix(A) do i,j,aij
        aij == 0 ? " " : nzmark
    end


"""
    print_lvector(v)

Print a LaTeX vector.
"""
print_lvector(v :: AbstractVector) =
    svector(v) do i,vi
        vi
    end


"""
    print_lmatrix(A)

Print a LaTeX matrix
"""
print_lmatrix(A :: AbstractMatrix) =
    smatrix(A) do i,j,aij
        aij
    end


"""
    c = black_cmatrix(m, n)

Create an m-by-n matrix of `"black"`
"""
function black_cmatrix(m, n)
    c = Matrix{String}(undef,m,n)
    c[:] .= "black"
    c
end


"""
    A, c = gray_zero_cmatrix!(A, c)

Rewrite entries corresponding to 0 in A with `"lightgray"` in c.
"""
function gray_zero_cmatrix!(A, c)
    c[A .== 0] .= "lightgray"
    A, c
end


"""
    c = gray_zero_cmatrix(A)

Produce a color matrix parallel to A where zero is `"lightgray"`,
nonzero is `"black"`.
"""
gray_zero_cmatrix(A) = gray_zero_cmatrix!(A, black_cmatrix(size(A)...))


"""
    print_cvector(v, [c])

Print a colored LaTeX vector (colors are in `c`).  If `c` is not provided,
default to black for nonzero and light gray for zero elements.
"""
print_cvector(v :: AbstractVector, c) =
    svector(v) do i,vi
        scolor(c[i], vi)
    end

print_cvector(v :: AbstractVector) =
    svector(v) do i,vi
        scolor(vi == 0 ? "lightgray" : "black", vi)
    end

"""
    print_cmatrix(A, c)

Print a colored LaTeX matrixr (colors are in `c`).  If `c` is not provided,
default to black for nonzero and light gray for zero elements.
"""
print_cmatrix(A :: AbstractMatrix, c) =
    smatrix(A) do i,j,aij
        scolor(c[i,j], aij)
    end

print_cmatrix(A :: AbstractMatrix) =
    smatrix(A) do i,j,aij
        scolor(aij == 0 ? "lightgray" : "black", aij)
    end

