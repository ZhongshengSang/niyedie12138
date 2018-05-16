points=(51)
timestep=(1)
crystal=(8) #(6 8 9)
diffusionAGG=(1)
diffusionXXX=(1)
c0_sat=(100.0) #(100.0 3.0)
DoD=(0.5 1.0 1.5 2.0 3.0 3.5 4.0)


for csat in "${c0_sat[@]}"; do
	# this loop changes the saturation concentration
	mkdir $csat

	for cr in "${crystal[@]}"; do
		# this for loop changes which crystal size we are using, 8 and 32
		mkdir $csat/$cr
	
		for ((i=1;i<2 ;i+=1)); do
			# this loop changes the electron equivalence we are looking at
			mkdir $csat/$cr/$i

				for nj in "${points[@]}";do
				# this loop changes the number of node points

					rm fortran_temp*

					sed s/ooo/$cr/g 		*_x.f95 			> 	fortran_temp_a.f95
					sed s/mmm/$i/g 			fortran_temp_a.f95  >   fortran_temp_b.f95
					sed s/dfdf/$diff/g 		fortran_temp_b.f95 	> 	fortran_temp_c.f95
					sed s/njnj/$nj/g 		fortran_temp_c.f95  >   fortran_temp_d.f95
					sed s/cmxx/$csat/g  	fortran_temp_d.f95  >	$csat/$cr/$i/Fe3O4_agg_Csat_Charge_6nm.f95
					cd $csat/$cr/$i/
					
					gfortran *.f95
					./a.out
					rm a.out
					rm *.mod
					
					cd ../../../
					sleep 0.1

				done
			
		done
	done

	rm  *temp*
done

