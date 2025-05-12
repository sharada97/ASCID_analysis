for n in 2 3 4
do
  time nohup python3 -u Simulation.py $n >> Figure/Out/DataTxt.txt
  echo "**********Done File $n**************"
done
