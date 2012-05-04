
NAME= report


all:
	latex $(NAME).tex
	#bibtex $(NAME)
	latex $(NAME).tex
	latex $(NAME).tex
	dvips $(NAME).dvi
	ps2pdf $(NAME).ps

clean:
	rm -f *.pdf *.out *.aux *.ps *.dvi *.log
