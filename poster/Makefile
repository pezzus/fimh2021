.PHONY=all

all: FIMH2021.pdf

FIMH2021.pdf: poster.pdf
	gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=$@ -dBATCH $^

poster.pdf: poster.tex
	latexmk -pdf $<

clean:
	latexmk -C poster
