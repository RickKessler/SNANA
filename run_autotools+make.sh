# run a few autotools commands to setup on new system or after git pull.

echo ' '
echo '# ============================================= '
echo '#     Begin autotools configure ... '
echo '# ============================================= '
echo ' '

autoreconf -fi     # -f=force -i = install
./configure

echo ' '
echo '# ============================================= '
echo '#     Begin make ... '
echo '# ============================================= '
echo ' '

cd src
make clean
make

