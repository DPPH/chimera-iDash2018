decrypt:
	../fhe/build/decrypt-program
	# results are in stat.txt 

cloud:
	rm -rf clouddir || echo
	mkdir clouddir
	#public/cloud keys from keygen phase
	cd clouddir; ln ../params.bin .
	cd clouddir; ln ../bk.key .
	cd clouddir; ln ../rk.key .
	cd clouddir; ln ../ks.key .
	cd clouddir; ln ../cloud_keyset.bin .
	#results from the encryption phase
	cd clouddir; ln ../dimensions.txt .
	cd clouddir; ln ../data.ctxt .
	cd clouddir; ln ../X_lvl2.bin .
	cd clouddir; ln ../S_lvl3.bin .
	cd clouddir; ln ../y_lvl2.bin .
	#run the cloud part
	cd clouddir; ../../log_regr/build/log_regr
	cd clouddir; ../../fhe/build/cloud-program
	#copy the result
	cp clouddir/numerator_lvl0.bin .
	cp clouddir/denominator_lvl0.bin .

cloud-inplace: 
	../log_regr/build/log_regr
	../fhe/build/cloud-program
	#../fhe/build/cloud-program-debug

encrypt: 
	../log_regr/build/encrypt_data
	../fhe/build/encrypt-program

keygen: 
	../log_regr/build/gen_keys
	../fhe/build/key-gen-program

compile: data compile_section1 compile_section2

compile_section1:
	bash ../log_regr/patch_tfhe_torus_generic.bash
	mkdir ../log_regr/build || echo
	cd ../log_regr/build; cmake .. -DCMAKE_BUILD_TYPE=Release
	cd ../log_regr/build; make -j 4

compile_section2:
	mkdir ../fhe/build || echo
	cd ../fhe/build; cmake .. -DCMAKE_BUILD_TYPE=Release
	cd ../fhe/build; make -j 4
