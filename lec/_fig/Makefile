FIGS=   subnormal0.svg \
	subnormal1.svg \
	error_fwd.svg \
	error_bwd.svg \
	gshear.svg \
	up-arrow-graph.svg \
	dn-arrow-graph.svg \
	nd-plot.svg

%.svg: out/%.xdv
	dvisvgm $<

out/%.xdv: %.tex
	xelatex --no-pdf --output-directory=out $<

all: $(FIGS)

subnormal0.svg subnormal1.svg: subnormal.jl
	julia subnormal.jl
