.PHONY: mex f2c

mex:
	/usr/local/MATLAB/R2019a/bin/mex gcvsplmex.c gcvspl.c
	/usr/local/MATLAB/R2019a/bin/mex spldermex.c gcvspl.c

f2c:
	@cp ./f2c/f2c.h .
	@./f2c/f2c gcvspl.f
	@echo "PATCH: replace comparison 'gf3 > gf2' by 'gf3 >= gf2' before comment 'Least-squares polynomial' to prevent stalls!"

