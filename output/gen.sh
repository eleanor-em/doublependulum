#!/bin/bash
echo '#!/bin/bash'                                                  >> run.sh
# resize the images
for i in `seq 0 499`;
do
    echo 'convert pgplot'$i'.png -resize 680x680\! pgplot'$i'.png'  >> run.sh
done
echo 'echo "images resized..."'                                     >> run.sh
# list all files in order to be added to the gif
echo -n 'convert -loop 0 '                                          >> run.sh
for i in `seq 0 499`;
do
    echo -n pgplot$i.png' '                                         >> run.sh
done
echo 'pendulum.gif'                                                 >> run.sh
echo 'echo "done!"'                                                 >> run.sh
echo "starting..."
./run.sh
# clean up temporary files
rm run.sh