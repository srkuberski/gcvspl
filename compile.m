clear( 'all' );

% alternative: use makefile

mex -v gcvsplmex.c gcvspl.c
mex -v spldermex.c gcvspl.c

