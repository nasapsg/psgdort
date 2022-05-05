all: psgtest

psgtest: psgdort.c psgtest.c
	gcc psgtest.c -o psgtest

clean:
	@rm -f core *a *o *BAK *bak *~ *% *.log psgtest
