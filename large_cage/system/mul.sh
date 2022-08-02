va=0
for z in 0 2 4 6 8
do
    for x in 0 5 10 15 20
    do
        cp sGx0z0 sGx$x-z$z 
        sed -i "s|(0|($x|g"  sGx$x-z$z  
        if [ $z != $va ]  
        then
            sed -i "s|-0.02)|-$z)|g"  sGx$x-z$z
        fi
        echo sGx$x-z$z 
    done
      
done


postProcess -funcs '(sGx0-z0 sGx5-z0 sGx10-z0 sGx15-z0 sGx20-z0 sGx0-z2 sGx5-z2 sGx10-z2 sGx15-z2 sGx20-z2 sGx0-z4 sGx5-z4 sGx10-z4 sGx15-z4 sGx20-z4 sGx0-z6 sGx5-z6 sGx10-z6 sGx15-z6 sGx20-z6 sGx0-z8 sGx5-z8 sGx10-z8 sGx15-z8 sGx20-z8)' -time '550:600'

