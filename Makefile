PANDOC = pandoc
TXT = exp/exp_doc.txt qr_decomp.txt quad_eqn_err.txt
OUT = $(TXT:.txt=.pdf)

FLAGS = \
  -V "mainfont:DejaVu Serif" \
  -V "monofont:DejaVu Sans Mono" \
  --pdf-engine=xelatex

all: $(OUT)

%.pdf: %.txt
	$(PANDOC) $(FLAGS) -f rst -t pdf -o $@ $<

# TODO: Support rst2pdf output nicely?
#report.pdf: gpu_report.rst
#	rst2pdf $^ -o $@

clean:
	rm -f $(OUT)
