rm time.txt

start=1
end=30

counter=1
while [ $counter -lt $end ]; do
#
sed -i "s/^nspectral =.*/nspectral =$counter/" main_matrix.py
#time="$(time ( python3 main_matrix.py ) 2>&1 1>/dev/null )"
utime="$( TIMEFORMAT='%lU';time ( python3 main_matrix.py ) 2>&1 1>/dev/null )"
cat >> time.txt <<EOL
${counter} ${utime}
EOL
#
echo "Done" ${counter} "Spectral modes, " ${end}-${counter} "to go"
let counter=counter+1
done
