all: render

render:
	quarto render --to all
	
html:
	quarto render