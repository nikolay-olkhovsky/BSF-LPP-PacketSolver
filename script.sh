for i in $(seq 0 9);
do
    mpiexec -n 8 .\\x64\\Debug\\BSF-LPP-PacketSolver.exe $i
done