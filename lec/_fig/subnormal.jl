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

text(f, p, attr, label) =
    println(f, "\\node at ($(p[1]),$(p[2])) [$attr] {$label};")

function plot_normalized(f)
    line(f, [-6,0], [6,0], "<->")
    line(f, [0,-1], [0,1])

    # 8-bit format: E5M2
    p = 3
    s = 1.5
    for binade = 0:1
        x = s*(binade+1)
        line(f, ( x,-0.4), ( x,0.4), "ultra thick")
        line(f, (-x,-0.4), (-x,0.4), "ultra thick")
        text(f, [s*(binade+1), -0.4], "below",
             "\\large \$2^{$(binade+1)-\\mathrm{bias}}\$")
        text(f, [-s*(binade+1), -0.4], "below",
             "\\large \$-2^{$(binade+1)-\\mathrm{bias}}\$")
        for k = 1:2^(p-1)-1
            x = s*((binade+1) + (k/2^(p-1))*(binade+1))
            line(f, ( x,-0.2), ( x,0.2), "thick")
            line(f, (-x,-0.2), (-x,0.2), "thick")
        end
    end
end

function plot_subnormals(f)
    p = 3
    s = 1.5
    for k = -2^(p-1)+1:2^(p-1)-1
        x = s*(k/2^(p-1))
        line(f, ( x,-0.2), ( x,0.2), "thick,red")
        line(f, (-x,-0.2), (-x,0.2), "thick,red")        
    end
end

open("subnormal0.tex", "w") do f
    header(f)
    plot_normalized(f)
    footer(f)
end

open("subnormal1.tex", "w") do f
    header(f)
    plot_normalized(f)
    plot_subnormals(f)
    footer(f)
end
