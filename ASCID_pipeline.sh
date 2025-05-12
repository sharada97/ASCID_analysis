for i in {0..2}
do
    for n in {0..3}
    do
      output_file="/mnt/e/RUN5/Take5/Ba133/Figures/Log/log_${i}_${n}.txt"
      time python3 -u ASCID_ThresholdEfficiency.py $i $n > "$output_file" 2>&1
      echo "**********Done File $n, folder $i**************"
      sleep 60
    done
done






# for i in 1 3 5 7 9
# do
#     for n in 6 8 10 12 14 16
#     do
#       output_file="../RUN5/Sc46/Figures/Take4_2/Log/log_${i}_${n}.txt"
#       time python3 -u ASCID_combine.py 9 1 $i $n > "$output_file" 2>&1
#       echo "**********Done File $n, folder $i**************"
#       sleep 60
#     done
# done





# i = 3
# n = 3
# output_file="../RUN5/Sc46/Figures/Take2/Log/log_${i}_${n}.txt"
# time nohup python3 -u ASCID_combine.py $i $n > "$output_file"

# output_file="../RUN5/Sc46/Figures/Take2/Log/Output_0_0.txt"
# time nohup python3 -u ASCID_fit.py 0 0 >> "$output_file"