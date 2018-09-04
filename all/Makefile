all:

cloud: X_beta.ctxt denominator_lvl0.bin

encrypt: data.ctxt S_lvl3.bin

keygen: cloud_keyset.bin bk.key

denominator_lvl0.bin: X_beta.ctxt
	../fhe/build/cloud-program

X_beta.ctxt:
	../log_regr/build/log_regr

S_lvl3.bin:
	../fhe/build/encrypt-program

data.ctxt:
	../log_regr/build/encrypt_data

bk.key: cloud_keyset.bin
	../fhe/build/key-gen-program


cloud_keyset.bin:
	../log_regr/build/gen_keys


compile: data compile_section2

compile_section2:
	mkdir ../fhe/build
	cd ../fhe/build; cmake .. -DCMAKE_BUILD_TYPE=Release
	cd ../fhe/build; make -j 4
