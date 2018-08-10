PATCH_DIR=$(cd `dirname $0` && pwd)
echo $PATCH_DIR
cd $PATCH_DIR/external/tfhe_torus_generic
git apply --check $PATCH_DIR/tfhe_torus_generic.patch && git am --signoff -k $PATCH_DIR/tfhe_torus_generic.patch
cd -
