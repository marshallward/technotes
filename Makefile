PANDOC = pandoc
OUT = out.pdf
TXT = quad_eqn_err.txt

FLAGS = \
  -V "mainfont:DejaVu Serif" \
  -V "monofont:DejaVu Sans Mono" \
  --pdf-engine=xelatex

all: $(OUT)

$(OUT): $(TXT)
	$(PANDOC) $(FLAGS) -f rst -t pdf -o $@ $<

# TODO: Support rst2pdf output nicely?
#report.pdf: gpu_report.rst
#	rst2pdf $^ -o $@

clean:
	rm -f $(OUT)
