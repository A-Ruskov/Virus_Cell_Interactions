for i in {1..1}
do
   mpirun -n 16  ./bin/Virus_Cell_Model.exe ./props/config.props ./props/model.props
   echo "Iteration $i completed!"
done