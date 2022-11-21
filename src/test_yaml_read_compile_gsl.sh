#export LGSL=$GSL_DIR/$LIB_SDIR/libgsl.a $GSL_DIR/$LIB_SDIR/libgslcblas.a
export LGSL=$GSL_DIR/lib/
export IGSL=-I$GSL_DIR/include
gcc -lyaml -lm -ldl -lgsl -lgslcblas $IGSL -L$LGSL -I/global/homes/g/gnarayan/local/include -L/global/homes/g/gnarayan/local/lib  test_yaml_read_gsl.c -o ./test_yaml_read_gsl
