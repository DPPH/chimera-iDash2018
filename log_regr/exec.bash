for iter in 0 1 2 3 4
do
	echo $iter
	./gen_keys
	./encrypt_data
	out_dir=out_$iter
	mkdir -p $out_dir
	/usr/bin/time -v --output=$out_dir/timing ./log_regr -vv --write_all_beta --write_all_x_beta --iters 10 > log_regr.out
	ls -l *ctxt *bin > $out_dir/sizes

	rm -f $out_dir/betas
	for f in beta*.ctxt; do echo $f `./decrypt_beta $f` >> $out_dir/betas ; done
	rm -f $out_dir/X_betas
	for f in X_beta*.ctxt; do echo $f `./decrypt_X_beta $f` >> $out_dir/X_betas ; done
done
