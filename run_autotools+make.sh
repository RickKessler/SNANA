# run a few autotools commands to setup on new system or after git pull.
# ./arg_configure        (normal configure+make)
# ./arg_configure debug  (debug  configure+make)
#

debug_flag=$1

if [[ $debug_flag == "debug" ]]; then
    arg_configure='--enable-debug'
else
    arg_configure=
fi


echo " "
echo "# ============================================= "
echo "#     Begin autotools configure $arg_configure "
echo "# ============================================= "
echo " "


autoreconf -fi     # -f=force -i = install
./configure $arg_configure

echo ' '
echo '# ============================================= '
echo '#     Begin make ... '
echo '# ============================================= '
echo ' '

cd src
make clean
make

