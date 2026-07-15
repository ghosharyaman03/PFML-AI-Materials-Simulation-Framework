cd Example_Systems/NiAlMo
mpirun --allow-run-as-root -np 4 /app/microsim_gp/microsim_gp Input_NiAlMo_4.in Filling_NiAlMo_4.in ./simulation_output 2 2
ls -lh ./simulation_output