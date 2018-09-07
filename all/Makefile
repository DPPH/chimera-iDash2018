decrypt:
	../fhe/build/decrypt-program
	# results are in stat.txt 

cloud:
	rm -rf cloud || echo
	mkdir cloud
	#public/cloud keys from keygen phase
	cd cloud; ln ../params.bin .
	cd cloud; ln ../bk.key .
	cd cloud; ln ../rk.key .
	cd cloud; ln ../ks.key .
	cd cloud; ln ../cloud_keyset.bin .
	#results from the encryption phase
	cd cloud; ln ../dimensions.txt .
	cd cloud; ln ../data.ctxt .
	cd cloud; ln ../X_lvl2.bin .
	cd cloud; ln ../S_lvl3.bin .
	cd cloud; ln ../y_lvl2.bin .
	#run the cloud part
	cd cloud; ../../log_regr/build/log_regr
	cd cloud; ../../fhe/build/cloud-program
	#copy the result
	cp cloud/numerator_lvl0.bin .
	cp cloud/denominator_lvl0.bin .

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
