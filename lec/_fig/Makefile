FIGS=   subnormal0.svg \
	subnormal1.svg \
	error_fwd.svg \
	error_bwd.svg \
	gshear.svg \
	up-arrow-graph.svg \
	up-arrow-graph-fill.svg \
	dn-arrow-graph.svg \
	nd-plot.svg \
	nd-plot0.svg \
	sgraph-A.svg \
	sgraph-A-1.svg \
	sgraph-A-2.svg \
	sgraph-A-3.svg \
	sgraph-A-4.svg \
	sgraph-A-5.svg \
	sgraph-A-6.svg \
	sgraph-A-7.svg \
	sgraph-A-8.svg \
	sgraph-A-9.svg \
	sgraph-A-10.svg \
	sgraph-A-11.svg \
	sgraph-A-12.svg \
	sgraph-A-13.svg \
	sgraph-A-14.svg \
	sgraph-L.svg \
	sgraph-LT.svg \
	lspic.svg \
	reflector.svg \
	qr-econ.svg \
	qr-QTr.svg \
	svd-econ.svg \
	svd-UTr.svg

%.svg: out/%.xdv
	dvisvgm -n $<

out/%.xdv: %.tex
	xelatex --no-pdf --output-directory=out $<

all: $(FIGS)

subnormal0.tex subnormal1.tex: subnormal.jl
	julia subnormal.jl

sgraph-A.tex sgraph-L.tex sgraph-LT.tex sgraph-AL-fig.tex: demo-sparse.jl
	julia demo-sparse.jl
