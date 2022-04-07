for i in {1..10}
do
   mpirun -n 64  ./bin/Virus_Cell_Model.exe ./props/config.props ./props/model.props
   echo "Iteration $i completed!"
done